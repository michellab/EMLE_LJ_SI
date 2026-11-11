"""Compute the B3LYP/cc-pVTZ wavefunction of each molecule with ORCA and write it as a molden file."""

import argparse
import os
import subprocess
from pathlib import Path

KEYWORDS = "B3LYP cc-pVTZ TightSCF DEFGRID3"


def main() -> None:
    """
    Run a single point for every xyz file in a directory.

    Returns
    -------
    None
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("geometries", type=Path, help="Directory with the xyz files (Angstrom)")
    parser.add_argument("output", type=Path, help="Output directory")
    parser.add_argument("--orca-dir", type=Path, required=True, help="ORCA installation directory")
    parser.add_argument("--nprocs", type=int, default=16)
    args = parser.parse_args()

    orca_dir, output = args.orca_dir.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PATH=f"{orca_dir}:{os.environ['PATH']}")

    for xyz_file in sorted(args.geometries.glob("*.xyz")):
        name = xyz_file.stem
        if (output / f"{name}.molden.input").exists():
            continue

        coordinates = "\n".join(xyz_file.read_text().splitlines()[2:]).strip()
        (output / f"{name}.inp").write_text(f"! {KEYWORDS}\n%pal nprocs {args.nprocs} end\n* xyz 0 1\n{coordinates}\n*\n")

        # ORCA must be called with its full path to run in parallel
        with open(output / f"{name}.out", "w") as f:
            subprocess.run([orca_dir / "orca", f"{name}.inp"], cwd=output, env=env, stdout=f, check=True)
        subprocess.run([orca_dir / "orca_2mkl", name, "-molden"], cwd=output, env=env, check=True, capture_output=True)
        print(f"{name} done", flush=True)


if __name__ == "__main__":
    main()
