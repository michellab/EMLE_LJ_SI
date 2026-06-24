"""Parse DES CSV files and generate .h5 files."""

import argparse
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from glob import glob
from pathlib import Path
from typing import Optional

import h5py
import numpy as np
import pandas as pd

try:
    import openmm as mm
    from openff.interchange import Interchange
    from openff.toolkit import ForceField, Molecule, Topology
    from openmm import unit
    from rdkit import Chem

    FORCE_FIELD_AVAILABLE = True
except ImportError:
    FORCE_FIELD_AVAILABLE = False


# Module-level globals for worker processes
_WORKER_DATA = {}


def _init_worker(molecules_data, xdm_data, force_field_file):
    """Initialize worker process with shared data."""
    global _WORKER_DATA
    _WORKER_DATA["xdm"] = xdm_data
    _WORKER_DATA["force_field"] = ForceField(force_field_file)
    _WORKER_DATA["molecules"] = {}
    _WORKER_DATA["rdkit_cache"] = {}
    _WORKER_DATA["charges_cache"] = {}

    for smiles, mol_data in molecules_data.items():
        rdkit_mol = Chem.MolFromMolBlock(mol_data, sanitize=False, removeHs=False)
        _WORKER_DATA["rdkit_cache"][smiles] = rdkit_mol
        _WORKER_DATA["molecules"][smiles] = Molecule.from_rdkit(
            rdkit_mol, allow_undefined_stereo=True, hydrogens_are_explicit=True
        )


def _get_charges_worker(solute_smiles, solvent_smiles):
    """Get MM charges in worker process."""
    key = (solute_smiles, solvent_smiles)
    if key not in _WORKER_DATA["charges_cache"]:
        topology = Topology.from_molecules(
            [
                _WORKER_DATA["molecules"][solute_smiles],
                _WORKER_DATA["molecules"][solvent_smiles],
            ]
        )
        interchange = Interchange.from_smirnoff(
            _WORKER_DATA["force_field"], topology, allow_nonintegral_charges=True
        )
        omm_system = interchange.to_openmm()
        nb = next(f for f in omm_system.getForces() if isinstance(f, mm.NonbondedForce))
        _WORKER_DATA["charges_cache"][key] = np.array(
            [
                nb.getParticleParameters(i)[0].value_in_unit(mm.unit.elementary_charge)
                for i in range(nb.getNumParticles())
            ]
        )
    return _WORKER_DATA["charges_cache"][key]


def _get_vdw_params_worker(solute_smiles, solvent_smiles):
    """Get vdW parameters (epsilon, sigma) for each atom in worker process."""
    key = (solute_smiles, solvent_smiles)
    if key not in _WORKER_DATA.get("vdw_cache", {}):
        if "vdw_cache" not in _WORKER_DATA:
            _WORKER_DATA["vdw_cache"] = {}

        topology = Topology.from_molecules(
            [
                _WORKER_DATA["molecules"][solute_smiles],
                _WORKER_DATA["molecules"][solvent_smiles],
            ]
        )
        interchange = Interchange.from_smirnoff(
            _WORKER_DATA["force_field"], topology, allow_nonintegral_charges=True
        )
        omm_system = interchange.to_openmm()
        nb = next(f for f in omm_system.getForces() if isinstance(f, mm.NonbondedForce))

        # Conversion factor: 1 nm = 18.897259886 Bohr
        NM_TO_BOHR = 18.897259886

        epsilons = []
        sigmas = []
        for i in range(nb.getNumParticles()):
            charge, sigma_nm, epsilon_kj = nb.getParticleParameters(i)
            # Convert sigma from nm to Bohr
            sigma_bohr = sigma_nm.value_in_unit(mm.unit.nanometer) * NM_TO_BOHR
            # Convert epsilon to kJ/mol
            epsilon_val = epsilon_kj.value_in_unit(mm.unit.kilojoule_per_mole)

            epsilons.append(epsilon_val)
            sigmas.append(sigma_bohr)

        _WORKER_DATA["vdw_cache"][key] = (np.array(epsilons), np.array(sigmas))

    return _WORKER_DATA["vdw_cache"][key]


def _process_system_worker(args):
    """Process a single system in worker process."""
    system_id, df_rows, solute_smiles, solvent_smiles = args

    KCAL_TO_KJ = 4.184

    solute_rdmol = _WORKER_DATA["rdkit_cache"][solute_smiles]
    solvent_rdmol = _WORKER_DATA["rdkit_cache"][solvent_smiles]

    solute_z = np.array(
        [a.GetAtomicNum() for a in solute_rdmol.GetAtoms()], dtype=np.int32
    )
    solvent_z = np.array(
        [a.GetAtomicNum() for a in solvent_rdmol.GetAtoms()], dtype=np.int32
    )
    n_solute = len(solute_z)

    charges = _get_charges_worker(solute_smiles, solvent_smiles)
    charges_mm = charges[n_solute:].astype(np.float64)

    # Get vdW parameters for all atoms
    epsilons_all, sigmas_all = _get_vdw_params_worker(solute_smiles, solvent_smiles)
    epsilon_qm = epsilons_all[:n_solute].astype(np.float64)
    sigma_qm = sigmas_all[:n_solute].astype(np.float64)
    epsilon_mm = epsilons_all[n_solute:].astype(np.float64)
    sigma_mm = sigmas_all[n_solute:].astype(np.float64)

    xdm_qm = _WORKER_DATA["xdm"].get(
        solute_smiles, (np.zeros(n_solute), np.zeros(n_solute), np.zeros(n_solute))
    )
    xdm_mm = _WORKER_DATA["xdm"].get(
        solvent_smiles,
        (np.zeros(len(solvent_z)), np.zeros(len(solvent_z)), np.zeros(len(solvent_z))),
    )

    expected = " ".join(
        [a.GetSymbol() for a in solute_rdmol.GetAtoms()]
        + [a.GetSymbol() for a in solvent_rdmol.GetAtoms()]
    )

    configs = []
    for row in df_rows:
        if row["elements"] != expected:
            continue

        coords = np.fromstring(row["xyz"], sep=" ", dtype=np.float64).reshape(-1, 3)

        # Convert all energy columns from kcal/mol to kJ/mol
        sapt_es = row["sapt_es"] * KCAL_TO_KJ
        sapt_ex = row["sapt_ex"] * KCAL_TO_KJ
        sapt_exs2 = row["sapt_exs2"] * KCAL_TO_KJ
        sapt_ind = row["sapt_ind"] * KCAL_TO_KJ
        sapt_exind = row["sapt_exind"] * KCAL_TO_KJ
        sapt_disp = row["sapt_disp"] * KCAL_TO_KJ
        sapt_exdisp_os = row["sapt_exdisp_os"] * KCAL_TO_KJ
        sapt_exdisp_ss = row["sapt_exdisp_ss"] * KCAL_TO_KJ
        sapt_delta_HF = row["sapt_delta_HF"] * KCAL_TO_KJ
        sapt_all = row["sapt_all"] * KCAL_TO_KJ
        cc_ccsd_t_all = row["cc_CCSD(T)_all"] * KCAL_TO_KJ

        configs.append(
            {
                "e_disp": sapt_disp + sapt_exdisp_os + sapt_exdisp_ss,
                "sapt_es": sapt_es,
                "sapt_ex": sapt_ex,
                "sapt_exs2": sapt_exs2,
                "sapt_ind": sapt_ind,
                "sapt_exind": sapt_exind,
                "sapt_disp": sapt_disp,
                "sapt_exdisp_os": sapt_exdisp_os,
                "sapt_exdisp_ss": sapt_exdisp_ss,
                "sapt_delta_HF": sapt_delta_HF,
                "sapt_all": sapt_all,
                "cc_ccsd_t_all": cc_ccsd_t_all,
                "xyz_qm": coords[:n_solute],
                "xyz_mm": coords[n_solute:],
                "xyz": coords,
                "atomic_numbers": solute_z,
                "atomic_numbers_mm": solvent_z,
                "charges_mm": charges_mm,
                "epsilon_qm": epsilon_qm,
                "sigma_qm": sigma_qm,
                "epsilon_mm": epsilon_mm,
                "sigma_mm": sigma_mm,
                "xdm_M1_sq": xdm_qm[0],
                "xdm_M2_sq": xdm_qm[1],
                "xdm_M3_sq": xdm_qm[2],
                "xdm_M1_sq_mm": xdm_mm[0],
                "xdm_M2_sq_mm": xdm_mm[1],
                "xdm_M3_sq_mm": xdm_mm[2],
            }
        )

    return system_id, configs


class DESCSVParser:
    """Parallelized parser for DES CSV files."""

    ALLOWED_ELEMENTS = frozenset(("H", "C", "N", "O", "S"))
    REQUIRED_COLUMNS = (
        "system_id",
        "smiles0",
        "smiles1",
        "charge0",
        "charge1",
        "natoms0",
        "natoms1",
        "xyz",
        "elements",
        "sapt_es",
        "sapt_ex",
        "sapt_exs2",
        "sapt_ind",
        "sapt_exind",
        "sapt_disp",
        "sapt_exdisp_os",
        "sapt_exdisp_ss",
        "sapt_delta_HF",
        "sapt_all",
        "cc_CCSD(T)_all",
    )

    def __init__(
        self,
        sdf_dir: str,
        force_field_file: str = "openff-2.0.0.offxml",
        xdm_data_dir: Optional[str] = None,
        n_workers: int = None,
        duplicate: bool = False,
    ):
        if not FORCE_FIELD_AVAILABLE:
            raise ImportError("Requires: rdkit, openff-toolkit, openmm")

        self.sdf_dir = Path(sdf_dir)
        self.force_field_file = force_field_file
        self.n_workers = n_workers or max(1, mp.cpu_count() - 1)
        self.duplicate = duplicate
        self.xdm_data = {}
        self.molecules_data = {}  # Store as mol blocks for pickling
        self._configs = []

        if xdm_data_dir:
            self._load_xdm_data(Path(xdm_data_dir))
        self._load_sdf_files()

    def _load_xdm_data(self, xdm_data_dir: Path) -> None:
        """Load XDM data from HDF5 files."""
        if not xdm_data_dir.exists():
            print(f"Warning: XDM data directory not found: {xdm_data_dir}")
            return

        xdm_files = glob(str(xdm_data_dir / "*/orca.molden.input.h5"))
        if not xdm_files:
            print(f"Warning: No XDM data files found in {xdm_data_dir}")
            return

        print(f"Loading XDM data from {len(xdm_files)} files")
        for fpath in xdm_files:
            try:
                smiles = Path(fpath).parent.name
                with h5py.File(fpath, "r") as f:
                    xdm = f["mbis"]["xdm"]["mbis"]["xdm_results"]
                    self.xdm_data[smiles] = (
                        xdm["<M1^2>"][:],
                        xdm["<M2^2>"][:],
                        xdm["<M3^2>"][:],
                    )
            except Exception as e:
                print(f"Warning: Could not load XDM from {fpath}: {e}")
        print(f"Loaded XDM for {len(self.xdm_data)} molecules")

    def _load_sdf_files(self) -> None:
        """Load SDF files as mol blocks (picklable)."""
        if not self.sdf_dir.exists():
            raise FileNotFoundError(f"SDF directory not found: {self.sdf_dir}")

        sdf_files = list(self.sdf_dir.glob("*.sdf"))
        if not sdf_files:
            raise ValueError(f"No SDF files found in {self.sdf_dir}")

        print(f"Loading {len(sdf_files)} SDF files")
        for sdf_file in sdf_files:
            try:
                smiles = sdf_file.stem
                rdkit_mol = Chem.SDMolSupplier(
                    str(sdf_file), sanitize=False, removeHs=False
                )[0]
                if rdkit_mol is None:
                    continue
                self.molecules_data[smiles] = Chem.MolToMolBlock(rdkit_mol)
            except Exception as e:
                print(f"Error loading {sdf_file}: {e}")

        print(f"Loaded {len(self.molecules_data)} molecules")
        if not self.molecules_data:
            raise ValueError("No valid molecules loaded")

    def process_csv(self, csv_file: str) -> None:
        """Process CSV file with parallel workers."""
        print(f"Processing {csv_file} with {self.n_workers} workers")
        df = pd.read_csv(csv_file)

        missing = set(self.REQUIRED_COLUMNS) - set(df.columns)
        if missing:
            raise ValueError(f"Missing columns: {missing}")

        # Tracking
        tasks = []
        skipped = {
            "unsupported_elements": [],
            "charged": [],
            "missing_sdf": [],
            "error": [],
        }
        kept = []

        grouped = df.groupby("system_id", sort=False)
        n_systems = len(grouped)
        print(f"Found {n_systems} systems in {len(df)} rows")

        for system_id, df_system in grouped:
            row = df_system.iloc[0]
            solute_smiles, solvent_smiles = row.smiles0, row.smiles1
            pair = f"{solute_smiles} + {solvent_smiles}"

            elements = set(row.elements.split())
            if not elements.issubset(self.ALLOWED_ELEMENTS):
                bad = elements - self.ALLOWED_ELEMENTS
                skipped["unsupported_elements"].append(
                    (system_id, pair, f"elements: {bad}")
                )
                continue

            if row.charge0 != 0 or row.charge1 != 0:
                skipped["charged"].append(
                    (system_id, pair, f"charges: {row.charge0}/{row.charge1}")
                )
                continue

            if solute_smiles not in self.molecules_data:
                skipped["missing_sdf"].append(
                    (system_id, pair, f"missing SDF: {solute_smiles}")
                )
                continue
            if solvent_smiles not in self.molecules_data:
                skipped["missing_sdf"].append(
                    (system_id, pair, f"missing SDF: {solvent_smiles}")
                )
                continue

            rows_data = df_system[
                [
                    "elements",
                    "xyz",
                    "sapt_es",
                    "sapt_ex",
                    "sapt_exs2",
                    "sapt_ind",
                    "sapt_exind",
                    "sapt_disp",
                    "sapt_exdisp_os",
                    "sapt_exdisp_ss",
                    "sapt_delta_HF",
                    "sapt_all",
                    "cc_CCSD(T)_all",
                ]
            ].to_dict("records")
            tasks.append((system_id, rows_data, solute_smiles, solvent_smiles))

        print(f"Submitting {len(tasks)} valid systems to workers...")

        valid_count = 0
        total_configs = 0

        with ProcessPoolExecutor(
            max_workers=self.n_workers,
            initializer=_init_worker,
            initargs=(self.molecules_data, self.xdm_data, self.force_field_file),
        ) as executor:
            futures = {
                executor.submit(_process_system_worker, task): task for task in tasks
            }

            for i, future in enumerate(as_completed(futures)):
                if i % 10 == 0 or i == len(futures) - 1:
                    print(
                        f"\rProcessed {i + 1}/{len(futures)} systems...",
                        end="",
                        flush=True,
                    )

                task = futures[future]
                system_id, _, solute, solvent = task
                pair = f"{solute} + {solvent}"

                try:
                    system_id, configs = future.result()
                    if configs:
                        self._configs.extend(configs)

                        # If duplicate flag is set, create swapped versions
                        if self.duplicate:
                            swapped_configs = self._create_swapped_configs(configs)
                            self._configs.extend(swapped_configs)

                        valid_count += 1
                        total_configs += len(configs) * (2 if self.duplicate else 1)
                        kept.append((system_id, pair, len(configs)))
                except Exception as e:
                    skipped["error"].append((system_id, pair, str(e)[:100]))
                    print(f"\n  Error: {system_id} ({pair})")

        # Print summary
        print(f"\n\n{'=' * 80}")
        print("PROCESSING SUMMARY")
        print(f"{'=' * 80}")
        print(f"\nKept: {valid_count} systems, {total_configs} configurations")

        total_skipped = sum(len(v) for v in skipped.values())
        print(f"Skipped: {total_skipped} systems")

        if skipped["unsupported_elements"]:
            print(f"\n  Unsupported elements ({len(skipped['unsupported_elements'])}):")
            for sid, pair, reason in skipped["unsupported_elements"][:10]:
                print(f"    {sid}: {pair} - {reason}")
            if len(skipped["unsupported_elements"]) > 10:
                print(f"    ... and {len(skipped['unsupported_elements']) - 10} more")

        if skipped["charged"]:
            print(f"\n  Charged molecules ({len(skipped['charged'])}):")
            for sid, pair, reason in skipped["charged"][:10]:
                print(f"    {sid}: {pair} - {reason}")
            if len(skipped["charged"]) > 10:
                print(f"    ... and {len(skipped['charged']) - 10} more")

        if skipped["missing_sdf"]:
            print(f"\n  Missing SDF ({len(skipped['missing_sdf'])}):")
            for sid, pair, reason in skipped["missing_sdf"][:10]:
                print(f"    {sid}: {pair} - {reason}")
            if len(skipped["missing_sdf"]) > 10:
                print(f"    ... and {len(skipped['missing_sdf']) - 10} more")

        if skipped["error"]:
            print(f"\n  Processing errors ({len(skipped['error'])}):")
            for sid, pair, reason in skipped["error"]:
                print(f"    {sid}: {pair}")
                print(f"      {reason}")

        print(f"\n{'=' * 80}")

    def _create_swapped_configs(self, configs):
        """Create duplicate configs with MM and QM regions swapped."""
        swapped = []
        for c in configs:
            swapped.append(
                {
                    # Energy values remain the same
                    "e_disp": c["e_disp"],
                    "sapt_es": c["sapt_es"],
                    "sapt_ex": c["sapt_ex"],
                    "sapt_exs2": c["sapt_exs2"],
                    "sapt_ind": c["sapt_ind"],
                    "sapt_exind": c["sapt_exind"],
                    "sapt_disp": c["sapt_disp"],
                    "sapt_exdisp_os": c["sapt_exdisp_os"],
                    "sapt_exdisp_ss": c["sapt_exdisp_ss"],
                    "sapt_delta_HF": c["sapt_delta_HF"],
                    "sapt_all": c["sapt_all"],
                    "cc_ccsd_t_all": c["cc_ccsd_t_all"],
                    # Swap QM and MM coordinates
                    "xyz_qm": c["xyz_mm"],
                    "xyz_mm": c["xyz_qm"],
                    "xyz": np.vstack([c["xyz_mm"], c["xyz_qm"]]),  # MM first, then QM
                    # Swap atomic numbers
                    "atomic_numbers": c["atomic_numbers_mm"],
                    "atomic_numbers_mm": c["atomic_numbers"],
                    # Swap charges (QM gets MM charges, MM gets zeros since original QM had no charges stored)
                    "charges_mm": np.zeros_like(c["atomic_numbers"], dtype=np.float64),
                    # Swap epsilon and sigma
                    "epsilon_qm": c["epsilon_mm"],
                    "sigma_qm": c["sigma_mm"],
                    "epsilon_mm": c["epsilon_qm"],
                    "sigma_mm": c["sigma_qm"],
                    # Swap XDM parameters
                    "xdm_M1_sq": c["xdm_M1_sq_mm"],
                    "xdm_M2_sq": c["xdm_M2_sq_mm"],
                    "xdm_M3_sq": c["xdm_M3_sq_mm"],
                    "xdm_M1_sq_mm": c["xdm_M1_sq"],
                    "xdm_M2_sq_mm": c["xdm_M2_sq"],
                    "xdm_M3_sq_mm": c["xdm_M3_sq"],
                }
            )
        return swapped

    def save_h5(self, output_file: str) -> None:
        """Save to HDF5 with optimized array stacking."""
        if not self._configs:
            raise ValueError("No data to save")

        n = len(self._configs)
        max_total = max(len(c["xyz"]) for c in self._configs)
        max_qm = max(len(c["atomic_numbers"]) for c in self._configs)
        max_mm = max(len(c["charges_mm"]) for c in self._configs)

        print(f"\nSaving {n} configs to {output_file}")
        print(f"Dims: total={max_total}, QM={max_qm}, MM={max_mm}")

        arrays = {
            "e_disp": np.empty(n, dtype=np.float64),
            "sapt_es": np.empty(n, dtype=np.float64),
            "sapt_ex": np.empty(n, dtype=np.float64),
            "sapt_exs2": np.empty(n, dtype=np.float64),
            "sapt_ind": np.empty(n, dtype=np.float64),
            "sapt_exind": np.empty(n, dtype=np.float64),
            "sapt_disp": np.empty(n, dtype=np.float64),
            "sapt_exdisp_os": np.empty(n, dtype=np.float64),
            "sapt_exdisp_ss": np.empty(n, dtype=np.float64),
            "sapt_delta_HF": np.empty(n, dtype=np.float64),
            "sapt_all": np.empty(n, dtype=np.float64),
            "cc_ccsd_t_all": np.empty(n, dtype=np.float64),
            "xyz": np.zeros((n, max_total, 3), dtype=np.float64),
            "xyz_qm": np.zeros((n, max_qm, 3), dtype=np.float64),
            "xyz_mm": np.zeros((n, max_mm, 3), dtype=np.float64),
            "atomic_numbers": np.zeros((n, max_qm), dtype=np.int32),
            "atomic_numbers_mm": np.zeros((n, max_mm), dtype=np.int32),
            "charges_mm": np.zeros((n, max_mm), dtype=np.float64),
            "epsilon_qm": np.zeros((n, max_qm), dtype=np.float64),
            "sigma_qm": np.zeros((n, max_qm), dtype=np.float64),
            "epsilon_mm": np.zeros((n, max_mm), dtype=np.float64),
            "sigma_mm": np.zeros((n, max_mm), dtype=np.float64),
            "xdm_M1_sq": np.zeros((n, max_qm), dtype=np.float64),
            "xdm_M2_sq": np.zeros((n, max_qm), dtype=np.float64),
            "xdm_M3_sq": np.zeros((n, max_qm), dtype=np.float64),
            "xdm_M1_sq_mm": np.zeros((n, max_mm), dtype=np.float64),
            "xdm_M2_sq_mm": np.zeros((n, max_mm), dtype=np.float64),
            "xdm_M3_sq_mm": np.zeros((n, max_mm), dtype=np.float64),
        }

        print("Building arrays...")
        for i, c in enumerate(self._configs):
            if i % 1000 == 0:
                print(f"\rBuilding arrays: {i}/{n} configs...", end="", flush=True)

            for key in (
                "e_disp",
                "sapt_es",
                "sapt_ex",
                "sapt_exs2",
                "sapt_ind",
                "sapt_exind",
                "sapt_disp",
                "sapt_exdisp_os",
                "sapt_exdisp_ss",
                "sapt_delta_HF",
                "sapt_all",
                "cc_ccsd_t_all",
            ):
                arrays[key][i] = c[key]

            ntot, nqm, nmm = len(c["xyz"]), len(c["xyz_qm"]), len(c["xyz_mm"])
            arrays["xyz"][i, :ntot] = c["xyz"]
            arrays["xyz_qm"][i, :nqm] = c["xyz_qm"]
            arrays["xyz_mm"][i, :nmm] = c["xyz_mm"]

            arrays["atomic_numbers"][i, :nqm] = c["atomic_numbers"]
            arrays["epsilon_qm"][i, :nqm] = c["epsilon_qm"]
            arrays["sigma_qm"][i, :nqm] = c["sigma_qm"]
            arrays["xdm_M1_sq"][i, :nqm] = c["xdm_M1_sq"]
            arrays["xdm_M2_sq"][i, :nqm] = c["xdm_M2_sq"]
            arrays["xdm_M3_sq"][i, :nqm] = c["xdm_M3_sq"]

            arrays["atomic_numbers_mm"][i, :nmm] = c["atomic_numbers_mm"]
            arrays["charges_mm"][i, :nmm] = c["charges_mm"]
            arrays["epsilon_mm"][i, :nmm] = c["epsilon_mm"]
            arrays["sigma_mm"][i, :nmm] = c["sigma_mm"]
            arrays["xdm_M1_sq_mm"][i, :nmm] = c["xdm_M1_sq_mm"]
            arrays["xdm_M2_sq_mm"][i, :nmm] = c["xdm_M2_sq_mm"]
            arrays["xdm_M3_sq_mm"][i, :nmm] = c["xdm_M3_sq_mm"]

        print(f"\rBuilding arrays: {n}/{n} configs... done")

        with h5py.File(output_file, "w") as f:
            # Save main datasets
            for key, arr in arrays.items():
                f.create_dataset(key, data=arr, compression="gzip", compression_opts=4)

            # Add attributes to document units for vdW parameters
            f["epsilon_qm"].attrs["units"] = "kJ/mol"
            f["sigma_qm"].attrs["units"] = "Bohr"
            f["epsilon_mm"].attrs["units"] = "kJ/mol"
            f["sigma_mm"].attrs["units"] = "Bohr"
            f.attrs["force_field"] = self.force_field_file

        print(
            f"Saved {len(arrays)} datasets (including per-atom OpenFF vdW parameters)"
        )


def main():
    parser = argparse.ArgumentParser(
        description="Parse DES CSV files for EMLE (parallelized)"
    )
    parser.add_argument("--csv", required=True, help="Input CSV file")
    parser.add_argument("--output", default="des_reference.h5", help="Output H5 file")
    parser.add_argument("--sdf-dir", required=True, help="SDF directory")
    parser.add_argument(
        "--force-field", default="openff-2.0.0.offxml", help="Force field"
    )
    parser.add_argument("--xdm-data-dir", default=None, help="XDM data directory")
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: CPU count - 1)",
    )
    parser.add_argument(
        "--duplicate",
        action="store_true",
        help="Duplicate data with MM and QM regions swapped",
    )

    args = parser.parse_args()

    p = DESCSVParser(
        sdf_dir=args.sdf_dir,
        force_field_file=args.force_field,
        xdm_data_dir=args.xdm_data_dir,
        n_workers=args.workers,
        duplicate=args.duplicate,
    )
    p.process_csv(args.csv)
    p.save_h5(args.output)


if __name__ == "__main__":
    main()
