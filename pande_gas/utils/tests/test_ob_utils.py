"""
Tests for ob_utils.
"""
import numpy as np
import unittest

from rdkit import Chem
from rdkit_utils import conformers

from pande_gas.utils import ob_utils


class TestIonizer(unittest.TestCase):
    """
    Test Ionizer.
    """
    def setUp(self):
        """
        Set up tests.
        """
        smiles = 'CC(C)Cc1ccc([C@H](C)C(=O)O)cc1'
        mol = Chem.MolFromSmiles(smiles)

        # generate conformers
        engine = conformers.ConformerGenerator(max_conformers=3)
        self.mol = engine.generate_conformers(mol)

        # ionize the oxygen manually
        ionized_mol = Chem.Mol(self.mol)  # create a copy
        found = False
        for atom in ionized_mol.GetAtoms():
            if atom.GetAtomicNum() == 8 and atom.GetTotalDegree() == 2:
                atom.SetFormalCharge(-1)
                atom.SetNoImplicit(True)
                found = True
        assert found
        self.ionized_mol = Chem.RemoveHs(ionized_mol)

        self.ionizer = ob_utils.Ionizer()

    def test_ionizer_flat(self):
        """
        Test Ionizer on molecules without 3D coordinates.
        """
        self.mol.RemoveAllConformers()
        assert self.mol.GetNumConformers() == 0

        # make sure Ionizer calls the right method
        try:
            self.ionizer._ionize_3d(self.mol)
            raise ValueError('Molecule should not have conformers.')
        except AssertionError:
            pass

        # compare molecule SMILES
        ionized_mol = self.ionizer(self.mol)
        assert Chem.MolToSmiles(
            ionized_mol, isomericSmiles=True) == Chem.MolToSmiles(
                self.ionized_mol, isomericSmiles=True)

    def test_ionizer_conformers(self):
        """
        Make sure ionization preserves heavy atom coordinates.
        """

        # make sure Ionizer calls the right method
        assert self.ionizer(self.mol).ToBinary() != self.ionizer._ionize_2d(
            self.mol).ToBinary()
        assert self.ionizer(self.mol).ToBinary() == self.ionizer._ionize_3d(
            self.mol).ToBinary()

        # compare molecule SMILES
        ionized_mol = self.ionizer(self.mol)
        ionized_mol = Chem.RemoveHs(ionized_mol)
        assert Chem.MolToSmiles(
            ionized_mol, isomericSmiles=True) == Chem.MolToSmiles(
                self.ionized_mol, isomericSmiles=True)

        # compare heavy atom coordinates
        assert ionized_mol.GetNumConformers() > 1  # multiple conformers
        assert (ionized_mol.GetNumConformers() ==
                self.ionized_mol.GetNumConformers())
        for a, b in zip(ionized_mol.GetConformers(),
                        self.ionized_mol.GetConformers()):
            for atom in ionized_mol.GetAtoms():
                idx = atom.GetIdx()
                a_pos = list(a.GetAtomPosition(idx))
                b_pos = list(b.GetAtomPosition(idx))
                b_pos = np.around(b_pos, 4)  # obabel rounds to four digits
                assert np.array_equal(a_pos, b_pos)

    def test_recombine_conformers(self):
        """
        Make sure Ionizer returns as many conformers as we give it.

        This enforces the assumption that a multiconformer molecule remains
        a multiconformer molecule even after ionization.
        """
        ionized_mol = self.ionizer(self.mol)
        assert ionized_mol.GetNumConformers() == self.mol.GetNumConformers()
        

class TestMolImage(unittest.TestCase):
    """
    Test MolImage.
    """
    def setUp(self):
        smiles = 'CC(C)Cc1ccc([C@H](C)C(=O)O)cc1'
        self.mol = Chem.MolFromSmiles(smiles)

    def test_images(self):
        """Test MolImage."""
        imager = ob_utils.MolImage()
        im = imager(self.mol)
        assert im.size == (32, 32), im.size
