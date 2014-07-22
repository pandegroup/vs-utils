"""
Test featurizer class.
"""
import numpy as np
import subprocess
import time
from unittest import skip

from rdkit import Chem
from rdkit_utils import conformers

from pande_gas.features.basic import MolecularWeight


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
    try:
        from IPython.parallel import Client
    except ImportError:
        skip('Cannot import from IPython.parallel.')
    cluster = subprocess.Popen(['ipcontroller', '--ip=*'])
    time.sleep(1)
    engine = subprocess.Popen(['ipengine', 'localhost'])
    time.sleep(10)
    mol = Chem.MolFromSmiles(test_smiles)
    f = MolecularWeight()
    rval = f([mol])
    parallel_rval = f([mol], parallel=True)
    assert np.array_equal(rval, parallel_rval)
    engine.kill()
    cluster.kill()

test_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
