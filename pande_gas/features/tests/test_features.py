"""
Test featurizer class.
"""
from rdkit import Chem

from pande_gas.features.basic import MolecularWeight
from pande_gas.utils import rdkit_utils as rd


def test_featurizer():
    """Test basic functionality of Featurizer."""
    mol = Chem.MolFromSmiles(test_smiles)
    f = MolecularWeight()
    rval = f([mol])
    assert rval.shape == (1, 1)


def test_flatten_conformers():
    """Flatten a multiconformer molecule."""
    mol = Chem.MolFromSmiles(test_smiles)
    mol = rd.generate_conformers(mol, 1)
    assert mol.GetNumConformers() > 0
    f = MolecularWeight()
    rval = f([mol])
    assert rval.shape == (1, 1)

test_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'

