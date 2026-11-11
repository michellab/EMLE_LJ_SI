"""Compare the Hirshfeld XDM moments, volumes and charges of PyXDM and postg.

For each ORCA molden file in a directory:

    postg                postg defaults; also writes its mesh and Hirshfeld weights (postg.patch)
    pyxdm_postg_inputs   PyXDM on postg's mesh and Hirshfeld weights
    pyxdm_postg_cutoff   the same, with b = 0 where the spin density is below 1e-10, as in postg
    pyxdm_<grid>         PyXDM defaults: HORTON grid and Hirshfeld weights from a pro-atom database

The results are written to a CSV file with one row per molecule, atom, quantity and method.
"""

import argparse
import csv
import os
import subprocess
from pathlib import Path

import horton as ht
import numpy as np
from pyxdm import XDMCalculator, XDMSession

HERE = Path(__file__).parent

QUANTITIES = ["<M1^2>", "<M2^2>", "<M3^2>", "volume", "charge"]


class PostgCutoffCalculator(XDMCalculator):
    """XDMCalculator with postg's cutoff: b = 0 where the spin density is below 1e-10."""

    def _compute_density_properties(self, grid_points: np.ndarray, anisotropic: bool = True) -> dict:
        props = super()._compute_density_properties(grid_points, anisotropic)
        for spin in ("alpha", "beta"):
            props[f"b_{spin}"] = np.where(props[f"rho_{spin}"] > 1e-10, props[f"b_{spin}"], 0.0)
        return props


def read_postg_table(lines: list[str], title: str, columns: slice) -> np.ndarray:
    """
    Read one per-atom table from the postg output.

    Parameters
    ----------
    lines : list[str]
        postg output
    title : str
        Line before the table header
    columns : slice
        Columns to read

    Returns
    -------
    np.ndarray
        Table of shape (n_atoms, n_columns)
    """
    start = lines.index(title) + 2
    end = lines.index("#", start)
    return np.array([[float(v) for v in line.split()[columns]] for line in lines[start:end]])


def run_postg(postg: Path, molden: Path, work_dir: Path) -> dict:
    """
    Run postg with its defaults, writing its mesh and Hirshfeld weights to work_dir.

    Parameters
    ----------
    postg : Path
        Patched postg executable (cleanmolden.awk is in the same directory)
    molden : Path
        ORCA molden file
    work_dir : Path
        Directory for postg's files

    Returns
    -------
    dict
        Results by quantity
    """
    name = molden.name.replace(".molden.input", "")
    with open(work_dir / f"{name}.molden", "w") as f:
        subprocess.run(["awk", "-f", postg.parent / "cleanmolden.awk", molden], stdout=f, check=True)

    # postg truncates its arguments to 100 characters, so it runs in work_dir with the file name only
    env = dict(os.environ, POSTG_DUMP_MESH=f"{name}.mesh.txt", POSTG_DUMP_WEIGHTS=f"{name}.weights.txt")
    command = [postg, "BJ", "0.0", "0.0", f"{name}.molden", "0.0"]
    output = subprocess.run(command, cwd=work_dir, env=env, capture_output=True, text=True, check=True).stdout
    (work_dir / f"{name}.postg.out").write_text(output)

    lines = output.splitlines()
    moments = read_postg_table(lines, "moments and volumes ", slice(2, 6))
    charges = read_postg_table(lines, "hirshfeld charges ", slice(2, 3))
    return dict(zip(QUANTITIES, [*moments.T, charges[:, 0]]))


def run_pyxdm(molden: Path, scheme: str, mesh: Path | None = None, grid: str = "fine", postg_cutoff: bool = False, **partition_args) -> dict:
    """
    Compute moments, volumes and charges with PyXDM.

    Parameters
    ----------
    molden : Path
        ORCA molden file
    scheme : str
        Partitioning scheme: "hirshfeld", or "external" for weights from a file
    mesh : Path, optional
        Mesh file; a HORTON grid of quality `grid` is used if None
    grid : str, optional
        HORTON grid quality (default is "fine")
    postg_cutoff : bool, optional
        Use postg's density cutoff for b (default is False)
    **partition_args
        proatomdb or weights_file

    Returns
    -------
    dict
        Results by quantity
    """
    session = XDMSession(str(molden))
    session.load_molecule()
    session.setup_grid(str(mesh) if mesh else None, grid_definition=grid)
    session.calculator = PostgCutoffCalculator(session.mol) if postg_cutoff else XDMCalculator(session.mol)
    session.setup_partition_schemes([scheme], **partition_args)

    part, calculator = session.partitions[scheme], session.calculator
    results = calculator.calculate_xdm_moments(part, session.grid, order=[1, 2, 3])[scheme]["xdm_results"]
    results["volume"] = calculator.calculate_radial_moments(part, order=3)[scheme]["radial_moments"]["<r^3>"]
    results["charge"] = session.partition_schemes[scheme].get_charges(session.mol, calculator.dm_full)[scheme]["charges"]
    return results


def main() -> None:
    """
    Run the comparison for every molden file in a directory.

    Returns
    -------
    None
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("wfn_dir", type=Path, help="Directory with the ORCA molden files (name.molden.input)")
    parser.add_argument("output", type=Path, help="Output CSV file")
    parser.add_argument("--postg", type=Path, required=True, help="Patched postg executable")
    parser.add_argument("--proatomdb", type=Path, default=HERE / "proatomdb" / "atoms.h5")
    parser.add_argument("--grids", nargs="+", default=["fine", "veryfine", "ultrafine"])
    parser.add_argument("--work", type=Path, default=HERE / "work")
    args = parser.parse_args()

    postg_exe = args.postg.resolve()
    args.work.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["molecule", "atom", "element", "quantity", "method", "value"])

        for molden in sorted(args.wfn_dir.glob("*.molden.input")):
            name = molden.name.replace(".molden.input", "")
            mesh, weights = args.work / f"{name}.mesh.txt", args.work / f"{name}.weights.txt"

            methods = {"postg": run_postg(postg_exe, molden, args.work)}
            methods["pyxdm_postg_inputs"] = run_pyxdm(molden, "external", mesh=mesh, weights_file=str(weights))
            methods["pyxdm_postg_cutoff"] = run_pyxdm(molden, "external", mesh=mesh, weights_file=str(weights), postg_cutoff=True)
            for grid in args.grids:
                methods[f"pyxdm_{grid}"] = run_pyxdm(molden, "hirshfeld", grid=grid, proatomdb=str(args.proatomdb))

            symbols = [ht.periodic[z].symbol for z in ht.IOData.from_file(str(molden)).numbers]
            for method, results in methods.items():
                for quantity in QUANTITIES:
                    for atom, value in enumerate(results[quantity]):
                        writer.writerow([name, atom, symbols[atom], quantity, method, repr(float(value))])
            f.flush()
            print(f"{name} done", flush=True)


if __name__ == "__main__":
    main()
