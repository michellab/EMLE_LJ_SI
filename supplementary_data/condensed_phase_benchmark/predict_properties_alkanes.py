from pathlib import Path
from typing import Optional

import descent.targets.thermo
import descent.utils.molecule
import numpy as np
import openmm
import smee
import smee.converters
import smee.mm
import torch
from descent.targets.thermo import SimulationConfig
from emle.models import EMLE
from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule, Topology
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper
from openff.units import unit as off_unit

MOLECULES = {
    "propane": {"smiles": "CCC", "temperature": 231.0},
    "butane": {"smiles": "CCCC", "temperature": 273.0},
    "hexane": {"smiles": "CCCCCC", "temperature": 261.0},
    "isobutane": {"smiles": "CC(C)C", "temperature": 298.0},
}
N_REPEATS = 3
NMOL = 1000
PRESSURE = 1.01325  # bar (1 atm)

ANGSTROM_TO_BOHR = 1.8897259886
HARTREE_TO_KJ_MOL = 2625.5002
EMLE_MODEL = (
    "../../training/02_fit_xdm_moments_reference/ligand_patched_static_only_xdm.mat"
)
BASE_FORCEFIELD = "openff-2.0.0.offxml"


def create_system_from_smiles(
    smiles_list: list[str],
    nmol_list: list[int],
    forcefield_name: str = "openff-2.0.0.offxml",
    charge_assignment_callback: Optional[callable] = None,
) -> tuple[smee.TensorSystem, smee.TensorForceField, list[smee.TensorTopology]]:
    """Create a tensor system from SMILES strings (one component per SMILES/nmol pair)."""
    force_field = ForceField(forcefield_name, load_plugins=True)

    try:
        mols = [
            Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)
            for smiles in smiles_list
        ]
    except Exception:
        mols = [
            Molecule.from_smiles(smiles, allow_undefined_stereo=True)
            for smiles in smiles_list
        ]

    if charge_assignment_callback is not None:
        for mol in mols:
            charge_assignment_callback(mol)

    interchanges = [
        Interchange.from_smirnoff(
            force_field,
            [mol],
            charge_from_molecules=[mol] if charge_assignment_callback else None,
        )
        for mol in mols
    ]
    tensor_forcefield, topologies = smee.converters.convert_interchange(interchanges)
    tensor_system = smee.TensorSystem(topologies, nmol_list, is_periodic=True)

    return tensor_system, tensor_forcefield, topologies


def build_simulation_config(temperature):
    """Build the bulk simulation config for a given temperature (K), always at 1 atm."""
    return {
        "bulk": SimulationConfig(
            max_mols=1000,
            gen_coords=smee.mm.GenerateCoordsConfig(),
            equilibrate=[
                smee.mm.SimulationConfig(
                    temperature=temperature * openmm.unit.kelvin,
                    pressure=None,
                    n_steps=50000,
                    timestep=1.0 * openmm.unit.femtosecond,
                ),
                smee.mm.SimulationConfig(
                    temperature=temperature * openmm.unit.kelvin,
                    pressure=PRESSURE * openmm.unit.bar,
                    n_steps=100000,
                    timestep=1.0 * openmm.unit.femtosecond,
                ),
            ],
            production=smee.mm.SimulationConfig(
                temperature=temperature * openmm.unit.kelvin,
                pressure=PRESSURE * openmm.unit.bar,
                n_steps=1000000,
                timestep=1.0 * openmm.unit.femtosecond,
            ),
            production_frequency=2000,
        )
    }


def build_emle_xdm_forcefield(smiles, device, output_path):
    """Patch openff-2.0.0's vdW types used by the molecule with EMLE/XDM LJ params fit to it."""
    mol = Molecule.from_smiles(smiles)
    mol.generate_conformers(n_conformers=1)
    topology = Topology.from_molecules([mol])

    xyz = torch.tensor(
        mol.conformers[0].m_as(off_unit.angstrom).astype(np.float32),
        dtype=torch.float32,
        device=device,
    ).unsqueeze(0)
    species = torch.tensor(
        [atom.atomic_number for atom in mol.atoms], dtype=torch.long, device=device
    ).unsqueeze(0)
    q_mol = torch.zeros(1, device=device)

    emle = EMLE(
        model=EMLE_MODEL, dispersion_mode="c6", device=device, alpha_mode="reference"
    )
    emle_base = emle._emle_base
    _, _, _, A_thole, c6, *_ = emle_base.forward(
        species, xyz, q_mol, calc_A_thole=True, calc_c6=True
    )
    alpha = emle_base.get_isotropic_polarizabilities_thole(A_thole).detach()
    c6 = c6.detach()

    c6 = 0.5 * c6 * alpha
    rmin = 2 * (2.54 * alpha ** (1 / 7))
    sigma = rmin / (2 ** (1 / 6))
    epsilon = c6 / (4 * sigma**6)
    sigma = (sigma / (ANGSTROM_TO_BOHR * 10)).cpu().numpy()[0]
    epsilon = (epsilon * HARTREE_TO_KJ_MOL).cpu().numpy()[0]

    force_field = ForceField(BASE_FORCEFIELD)
    vdw_handler = force_field.get_parameter_handler("vdW")
    labels = force_field.label_molecules(topology)[0]["vdW"]

    per_type_sigma, per_type_epsilon = {}, {}
    for atom_indices, param in labels.items():
        i = atom_indices[0]
        per_type_sigma.setdefault(param.id, []).append(sigma[i])
        per_type_epsilon.setdefault(param.id, []).append(epsilon[i])

    for param in vdw_handler.parameters:
        if param.id not in per_type_sigma:
            continue
        param.sigma = float(np.mean(per_type_sigma[param.id])) * off_unit.nanometer
        param.epsilon = (
            float(np.mean(per_type_epsilon[param.id])) * off_unit.kilojoules_per_mole
        )

    force_field.to_file(output_path)
    return output_path


def print_vdw_params(label, forcefield_path, topology):
    """Print sigma (Angstrom) / epsilon (kJ/mol) vdW params for the molecule's atom types."""
    force_field = ForceField(forcefield_path)
    labels = force_field.label_molecules(topology)[0]["vdW"]

    seen = {}
    for atom_indices, param in labels.items():
        seen.setdefault(param.id, param)

    print(f"{label} vdW parameters (molecule atom types):")
    for type_id in sorted(seen):
        param = seen[type_id]
        sigma_a = param.sigma.m_as(off_unit.angstrom)
        epsilon_kjmol = param.epsilon.m_as(off_unit.kilojoule_per_mole)
        print(
            f"  {type_id:<6} sigma = {sigma_a:8.4f} A   epsilon = {epsilon_kjmol:8.4f} kJ/mol"
        )


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

for name, molecule_info in MOLECULES.items():
    smiles = molecule_info["smiles"]
    temperature = molecule_info["temperature"]
    print(f"=== {name} ({smiles}) T = {temperature} K, P = 1 atm ===")

    simulation_config = build_simulation_config(temperature)

    emle_xdm_forcefield = build_emle_xdm_forcefield(
        smiles, device, f"openff-2.0.0-emle-xdm-{name}.offxml"
    )
    forcefields = {"OpenFF-2.0.0": BASE_FORCEFIELD, "EMLE/XDM": emle_xdm_forcefield}

    mapped_smiles = descent.utils.molecule.map_smiles(smiles)
    molecule_topology = Topology.from_molecules([Molecule.from_smiles(smiles)])

    # Reference value/std are unused for prediction, so are left as zero.
    dataset = descent.targets.thermo.create_dataset(
        {
            "type": "density",
            "smiles_a": smiles,
            "x_a": 1.0,
            "smiles_b": None,
            "x_b": None,
            "temperature": temperature,
            "pressure": PRESSURE,
            "value": 0.0,
            "std": 0.0,
            "units": "g/mL",
            "source": None,
        },
        {
            "type": "hvap",
            "smiles_a": smiles,
            "x_a": 1.0,
            "smiles_b": None,
            "x_b": None,
            "temperature": temperature,
            "pressure": PRESSURE,
            "value": 0.0,
            "std": 0.0,
            "units": "kcal/mol",
            "source": None,
        },
    )

    for label, forcefield in forcefields.items():
        print(f"Running {label} ({forcefield}) ...")
        print_vdw_params(label, forcefield, molecule_topology)

        densities, hvaps = [], []
        for repeat in range(N_REPEATS):
            _, tensor_ff, topologies = create_system_from_smiles(
                smiles_list=[smiles], nmol_list=[NMOL], forcefield_name=forcefield
            )
            tensor_ff = tensor_ff.to(device)
            topology_map = {mapped_smiles: topologies[0].to(device)}

            _, _, pred_vals, _ = descent.targets.thermo.predict(
                dataset,
                tensor_ff,
                topology_map,
                output_dir=Path(f"output/{name}/{label}/repeat_{repeat}/predictions"),
                cached_dir=Path(f"output/{name}/{label}/repeat_{repeat}/cache"),
                verbose=True,
                simulation_config=simulation_config,
            )
            densities.append(float(pred_vals[0]))
            hvaps.append(float(pred_vals[1]) * 4.184)  # kcal/mol -> kJ/mol

        density_mean, density_std = np.mean(densities), np.std(densities, ddof=1)
        hvap_mean, hvap_std = np.mean(hvaps), np.std(hvaps, ddof=1)
        print(
            f"{name} | {label}: density = {density_mean:.4f} +/- {density_std:.4f} g/mL "
            f"(n={N_REPEATS})"
        )
        print(
            f"{name} | {label}: Hvap    = {hvap_mean:.4f} +/- {hvap_std:.4f} kJ/mol "
            f"(n={N_REPEATS})"
        )
