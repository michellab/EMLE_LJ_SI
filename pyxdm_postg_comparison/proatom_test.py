"""PyXDM on a HORTON grid with Hirshfeld weights from postg's free-atom densities.

postg stores r^2 rho(r) of each free atom on the mesh q = r / (r + r_mid), r_mid = Z^(-1/3), and interpolates it with a
natural cubic spline in q (atomicdata.f90, and atomin in wfnmod.f90). The densities are read from postg's source and
interpolated the same way. The grid is the same as in the pyxdm_<grid> runs of compare.py, only the pro-atoms differ:
the difference from pyxdm_<grid> is the effect of the pro-atoms, and the difference from postg what is left.

The results are written to a CSV file in the format of compare.py, with the method pyxdm_<grid>_postg_proatoms.
"""

import argparse
import csv
import re
from pathlib import Path

import horton as ht
import numpy as np
from pyxdm import XDMCalculator, XDMSession
from pyxdm.partitioning.partitioning import ExternalPartition, ExternalPartitioning
from scipy.interpolate import CubicSpline

QUANTITIES = ["<M1^2>", "<M2^2>", "<M3^2>", "volume", "charge"]


def read_postg_densities(atomicdata: Path, elements: set[int]) -> dict[int, np.ndarray]:
    """
    Read the free-atom densities from postg's atomicdata.f90.

    Parameters
    ----------
    atomicdata : Path
        postg's atomicdata.f90
    elements : set[int]
        Atomic numbers to read

    Returns
    -------
    dict[int, np.ndarray]
        r^2 rho(r) on postg's mesh, points 1 to n, by atomic number
    """
    text = atomicdata.read_text()
    tables = {}
    for start, end, z, block in re.findall(r"ftot\(\s*(\d+):\s*(\d+),\s*(\d+)\)\s*=\s*\(/(.*?)/\)", text, flags=re.S):
        z = int(z)
        if z not in elements:
            continue
        values = [float(v) for v in block.replace("&", " ").replace(",", " ").split()]
        assert len(values) == int(end) - int(start) + 1
        tables.setdefault(z, {})[int(start)] = values
    return {z: np.concatenate([parts[k] for k in sorted(parts)]) for z, parts in tables.items()}


def free_atom_density(f: np.ndarray, z: int, r: np.ndarray) -> np.ndarray:
    """
    Free-atom density of postg at distances r, as in its subroutine atomin.

    Parameters
    ----------
    f : np.ndarray
        r^2 rho(r) on postg's mesh, points 1 to n
    z : int
        Atomic number
    r : np.ndarray
        Distances from the nucleus (Bohr)

    Returns
    -------
    np.ndarray
        rho(r)
    """
    n = len(f)
    # Knots q = j / (n + 1), j = 0 to n + 1, with r^2 rho = 0 at both ends (r = 0 and r = infinity)
    q_knots = np.arange(n + 2) / (n + 1)
    spline = CubicSpline(q_knots, np.concatenate([[0.0], f, [0.0]]), bc_type="natural")
    r_mid = z ** (-1.0 / 3.0)
    return np.abs(spline(r / (r + r_mid))) / r**2


def postg_weights(numbers: np.ndarray, coordinates: np.ndarray, points: np.ndarray, densities: dict[int, np.ndarray]) -> np.ndarray:
    """
    Hirshfeld weights from postg's free-atom densities.

    Parameters
    ----------
    numbers : np.ndarray
        Atomic numbers
    coordinates : np.ndarray
        Nuclear coordinates (Bohr), shape (n_atoms, 3)
    points : np.ndarray
        Points (Bohr), shape (n_points, 3)
    densities : dict[int, np.ndarray]
        From read_postg_densities

    Returns
    -------
    np.ndarray
        Weights, shape (n_points, n_atoms)
    """
    rho = np.empty((len(points), len(numbers)))
    for i, (z, centre) in enumerate(zip(numbers, coordinates)):
        rho[:, i] = free_atom_density(densities[z], z, np.linalg.norm(points - centre, axis=1))
    return rho / np.maximum(rho.sum(axis=1, keepdims=True), 1e-40)


def run(molden: Path, grid: str, densities: dict) -> dict:
    """
    Compute moments, volumes and charges with PyXDM on a HORTON grid and postg's pro-atoms.

    Parameters
    ----------
    molden : Path
        ORCA molden file
    grid : str
        HORTON grid quality
    densities : dict
        postg's free-atom densities, from read_postg_densities

    Returns
    -------
    dict
        Results by quantity
    """
    session = XDMSession(str(molden))
    session.load_molecule()
    session.setup_grid(None, grid_definition=grid)
    mol, molgrid = session.mol, session.grid
    calculator = XDMCalculator(mol)

    scheme = ExternalPartitioning()
    scheme._partition_obj = ExternalPartition(molgrid, postg_weights(mol.numbers, mol.coordinates, molgrid.points, densities))
    part = scheme.get_partition_object()

    results = calculator.calculate_xdm_moments(part, molgrid, order=[1, 2, 3])["external"]["xdm_results"]
    results["volume"] = calculator.calculate_radial_moments(part, order=3)["external"]["radial_moments"]["<r^3>"]
    results["charge"] = scheme.get_charges(mol, calculator.dm_full)["external"]["charges"]
    return results


def main() -> None:
    """
    Run every molden file given and write the CSV.

    Returns
    -------
    None
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("moldens", type=Path, nargs="+", help="ORCA molden files (name.molden.input)")
    parser.add_argument("--output", type=Path, required=True, help="Output CSV file")
    parser.add_argument("--atomicdata", type=Path, required=True, help="postg's atomicdata.f90")
    parser.add_argument("--grid", default="ultrafine", help="HORTON grid quality (default is ultrafine)")
    args = parser.parse_args()

    densities = read_postg_densities(args.atomicdata, {1, 6, 7, 8, 16})
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["molecule", "atom", "element", "quantity", "method", "value"])
        for molden in args.moldens:
            name = molden.name.replace(".molden.input", "")
            results = run(molden, args.grid, densities)
            symbols = [ht.periodic[z].symbol for z in ht.IOData.from_file(str(molden)).numbers]
            for quantity in QUANTITIES:
                for atom, value in enumerate(results[quantity]):
                    writer.writerow([name, atom, symbols[atom], quantity, f"pyxdm_{args.grid}_postg_proatoms", repr(float(value))])
            f.flush()
            print(f"{name} done", flush=True)


if __name__ == "__main__":
    main()
