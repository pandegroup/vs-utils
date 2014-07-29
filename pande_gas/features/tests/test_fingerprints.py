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

    def test_circular_fingerprints(self):
        """Test CircularFingerprint."""
        f = fp.CircularFingerprint()
        rval = f([self.mol])
        assert rval.shape == (1, f.size)
