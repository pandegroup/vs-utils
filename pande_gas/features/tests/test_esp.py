"""
Tests for electrostatic potential features.
"""
from rdkit import Chem
from rdkit_utils import conformers

from pande_gas.features.esp import ESP


def test_esp():
    """Test ESP."""
    mol = Chem.MolFromSmiles(test_smiles)
    mol = conformers.generate_conformers(mol)
    f = ESP()
    rval = f([mol])

    # the features array should contain:
    # * one molecule (first index)
    # * one conformer (second index)
    # * one cubic grid (remaining indices)
    assert rval.shape[:2] == (1, 1)
    size = rval.shape[2]
    assert rval.shape[2:] == (size, size, size)

test_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
