"""
Test topological fingerprints.
"""
import unittest

from rdkit import Chem

from pande_gas.features import fingerprints as fp


class TestCircularFingerprint(unittest.TestCase):
    """
    Tests for CircularFingerprint.
    """
    def setUp(self):
        """
        Set up tests.
        """
        smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
        self.mol = Chem.MolFromSmiles(smiles)
        self.engine = fp.CircularFingerprint()

    def test_circular_fingerprints(self):
        """
        Test CircularFingerprint.
        """
        rval = self.engine([self.mol])
        assert rval.shape == (1, self.engine.size)

    def test_sparse_circular_fingerprints(self):
        """
        Test CircularFingerprint with sparse encoding.
        """
        self.engine = fp.CircularFingerprint(sparse=True)
        rval = self.engine([self.mol])
        assert rval.shape == (1,)
        assert isinstance(rval[0], dict)
        assert len(rval[0])
