"""
Tests for electrostatic potential features.
"""
from rdkit import Chem
from rdkit_utils import conformers

from pande_gas.features.esp import ESP


def test_esp():
    """Test ESP."""
    mols = [Chem.MolFromSmiles(smiles) for smiles in test_smiles]
    mols = [conformers.generate_conformers(mol) for mol in mols]
    f = ESP()
    rval = f(mols)

    # the features array should contain:
    # * two molecules (first index)
    # * one conformer per molecule (second index)
    # * cubic grids for each conformer (remaining indices)
    assert rval.shape[:2] == (2, 1)
    size = rval.shape[2]
    assert rval.shape[2:] == (size, size, size)

test_smiles = ['CC(=O)OC1=CC=CC=C1C(=O)O', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O']
