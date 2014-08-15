"""
Test image featurizer.
"""
import unittest

from rdkit import Chem

from pande_gas.features import images


class TestOBabelMolImage(unittest.TestCase):
    """
    Test MolImage using the 'obabel' engine.
    """
    def setUp(self):
        """
        Set up tests.
        """
        smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
        self.mol = Chem.MolFromSmiles(smiles)
        self.engine = 'obabel'

    def test_images(self):
        """
        Test MolImage.
        """
        f = images.MolImage(250, flatten=True, engine=self.engine)
        rval = f([self.mol])
        assert rval.shape == (1, 250 * 250 * 3), rval.shape

    def test_images_topo_view(self):
        """
        Test MolImage using topo_view.
        """
        f = images.MolImage(250, flatten=False, engine=self.engine)
        rval = f([self.mol])
        assert rval.shape == (1, 250, 250, 3), rval.shape


class TestRDKitMolImage(TestOBabelMolImage):
    """
    Test MolImage using the 'rdkit' engine.
    """
    def setUp(self):
        """
        Set up tests.
        """
        super(TestRDKitMolImage, self).setUp()
        self.engine = 'rdkit'
