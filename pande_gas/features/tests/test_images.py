"""
Test image featurizer.
"""
import unittest

from rdkit import Chem

from pande_gas.features import images


class TestMolImage(unittest.TestCase):
    """
    Test MolImage.
    """
    def setUp(self):
        """
        Set up tests.
        """
        smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
        self.mol = Chem.MolFromSmiles(smiles)

    def test_images(self):
        """
        Test MolImage.
        """
        f = images.MolImage(250, flatten=True)
        rval = f([self.mol])
        assert rval.shape == (1, 250 * 250 * 3), rval.shape

    def test_images_topo_view(self):
        """
        Test MolImage using topo_view.
        """
        f = images.MolImage(250, flatten=False)
        rval = f([self.mol])
        assert rval.shape == (1, 250, 250, 3), rval.shape
