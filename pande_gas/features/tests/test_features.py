"""
Test featurizer class.
"""
import numpy as np

from rdkit import Chem
from rdkit_utils import conformers

from pande_gas.features.basic import MolecularWeight
from pande_gas.utils.parallel_utils import LocalCluster


def test_featurizer():
    """Test basic functionality of Featurizer."""
    mol = Chem.MolFromSmiles(test_smiles)
    f = MolecularWeight()
    rval = f([mol])
    assert rval.shape == (1, 1)


def test_flatten_conformers():
    """Flatten a multiconformer molecule."""
    mol = Chem.MolFromSmiles(test_smiles)
    mol = conformers.generate_conformers(mol, 1)
    assert mol.GetNumConformers() > 0
    f = MolecularWeight()
    rval = f([mol])
    assert rval.shape == (1, 1)


def test_parallel():
    """Test parallel featurization."""
    cluster = LocalCluster(1)
    mol = Chem.MolFromSmiles(test_smiles)
    f = MolecularWeight()
    rval = f([mol])
    parallel_rval = f([mol], parallel=True,
                      client_kwargs={'cluster_id': cluster.cluster_id})
    assert np.array_equal(rval, parallel_rval)

test_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
