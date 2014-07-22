"""
Test featurizer class.
"""
import numpy as np
import subprocess
import time
import uuid
from unittest import skip, TestCase

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


class TestParallel(TestCase):
    def setUp(self):
        try:
            from IPython.parallel import Client
        except ImportError:
            skip('Cannot import from IPython.parallel.')
        self.cluster_id = uuid.uuid4()
        self.controller = subprocess.Popen(
            ['ipcontroller', '--ip=*',
             '--cluster-id={}'.format(self.cluster_id)])
        time.sleep(1)
        self.engines = subprocess.Popen(
            ['ipengine', '--cluster-id={}'.format(self.cluster_id)])
        time.sleep(10)

    def test_parallel(self):
        """Test parallel featurization."""
        mol = Chem.MolFromSmiles(test_smiles)
        f = MolecularWeight()
        rval = f([mol])
        parallel_rval = f([mol], parallel=True,
                          client_kwargs={'cluster_id': self.cluster_id})
        assert np.array_equal(rval, parallel_rval)

    def tearDown(self):
        self.controller.terminate()
        self.engines.terminate()

test_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
