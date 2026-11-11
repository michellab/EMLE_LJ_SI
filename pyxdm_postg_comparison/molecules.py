"""Test molecules for the PyXDM-postg comparison and their starting geometries."""

import sys
from pathlib import Path

MOLECULES = {
    "n-pentane": "CCCCC",
    "n-octane": "CCCCCCCC",
    "benzene": "c1ccccc1",
    "naphthalene": "c1ccc2ccccc2c1",
    "water": "O",
    "methanol": "CO",
    "o-cresol": "Cc1ccccc1O",
    "aceticacid": "CC(=O)O",
    "1_4-dioxane": "C1COCCO1",
    "methylpentanoate": "CCCCC(=O)OC",
    "methylbenzoate": "COC(=O)c1ccccc1",
    "2-pentanone": "CCCC(C)=O",
    "ethylphenylether": "CCOc1ccccc1",
    "m-hydroxybenzaldehyde": "O=Cc1cccc(O)c1",
    "acetonitrile": "CC#N",
    "pyridine": "c1ccncc1",
    "aniline": "Nc1ccccc1",
    "benzamide": "NC(=O)c1ccccc1",
    "1-propanethiol": "CCCS",
    "thioanisole": "CSc1ccccc1",
}


def write_geometries(output_dir: Path) -> None:
    """
    Write a force-field geometry of each molecule as an xyz file (Angstrom).

    Parameters
    ----------
    output_dir : Path
        Directory for the xyz files

    Returns
    -------
    None
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    output_dir.mkdir(parents=True, exist_ok=True)
    for name, smiles in MOLECULES.items():
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol, randomSeed=7)
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol)
        else:
            AllChem.UFFOptimizeMolecule(mol)
        Chem.MolToXYZFile(mol, str(output_dir / f"{name}.xyz"))


if __name__ == "__main__":
    write_geometries(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "geometries")
