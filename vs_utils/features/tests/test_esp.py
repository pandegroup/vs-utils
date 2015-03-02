"""
Tests for electrostatic potential features.
"""
import unittest

from rdkit import Chem

from vs_utils.features.esp import ESP
from vs_utils.utils.rdkit_utils import conformers


class TestESP(unittest.TestCase):
    """
    Tests for ESP.
    """
    def setUp(self):
        """
        Set up tests.
        """
        smiles = ['CC(=O)OC1=CC=CC=C1C(=O)O', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O']
        mols = [Chem.MolFromSmiles(s) for s in smiles]
        engine = conformers.ConformerGenerator(max_conformers=1)
        self.mols = [engine.generate_conformers(mol) for mol in mols]
        for mol in self.mols:
            assert mol.GetNumConformers() > 0

    def test_esp(self):
        """
        Test ESP.
        """
        f = ESP()
        rval = f(self.mols)

        # the features array should contain:
        # * two molecules (first index)
        # * one conformer per molecule (second index)
        # * cubic grids for each conformer (remaining indices)
        assert rval.shape[:2] == (2, 1)
        size = rval.shape[2]
        assert rval.shape[2:] == (size, size, size)
