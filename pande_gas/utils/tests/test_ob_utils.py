"""
Tests for ob_utils.
"""
from rdkit import Chem

from pande_gas.utils import ob_utils


def test_ionizer():
    """Test Ionizer."""
    mol = Chem.MolFromSmiles(test_smiles)
    ionizer = ob_utils.Ionizer()
    ref_smiles = 'CC(C)Cc1ccc(C(C)C(=O)[O-])cc1'
    ionized_mol = ionizer(mol)
    ionized_smiles = Chem.MolToSmiles(ionized_mol, isomericSmiles=True,
                                      canonical=True)
    assert ionized_smiles == ref_smiles, ionized_smiles


def test_images():
    """Test MolImage."""
    mol = Chem.MolFromSmiles(test_smiles)
    imager = ob_utils.MolImage()
    im = imager(mol)
    assert im.size == (32, 32), im.size

test_smiles = 'CC(C)Cc1ccc(C(C)C(=O)O)cc1'
