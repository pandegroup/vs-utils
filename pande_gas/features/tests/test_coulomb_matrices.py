"""
Tests for Coulomb matrix calculation.
"""
import numpy as np

from rdkit import Chem
from rdkit_utils import conformers

from pande_gas.features import coulomb_matrices as cm


def test_coulomb_matrix():
    """Test CoulombMatrix."""
    mol = Chem.MolFromSmiles(test_smiles)
    mol = conformers.generate_conformers(mol, n_conformers=1)
    f = cm.CoulombMatrix(mol.GetNumAtoms())
    rval = f([mol])
    size = np.triu_indices(mol.GetNumAtoms())[0].size
    assert rval.shape == (1, 1, size)


def test_coulomb_matrix_padding():
    """Test CoulombMatrix with padding."""
    mol = Chem.MolFromSmiles(test_smiles)
    mol = conformers.generate_conformers(mol, n_conformers=1)
    f = cm.CoulombMatrix(max_atoms=mol.GetNumAtoms()*2)
    rval = f([mol])
    size = np.triu_indices(mol.GetNumAtoms()*2)[0].size
    assert rval.shape == (1, 1, size)

test_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
