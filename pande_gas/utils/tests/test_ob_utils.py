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
        engine = conformers.ConformerGenerator(max_conformers=1)
        smiles = 'CC(C)Cc1ccc([C@H](C)C(=O)O)cc1'
        mol = Chem.MolFromSmiles(smiles)
        self.mol = engine.generate_conformers(mol)

        # ionize the oxygen manually
        self.ionized_mol = Chem.Mol(self.mol)  # create a copy
        found = False
        for atom in self.ionized_mol.GetAtoms():
            if atom.GetAtomicNum() == 8 and atom.GetDegree() == 2:
                atom.SetFormalCharge(-1)
                found = True
        assert found

        self.ionizer = ob_utils.Ionizer()

    def test_ionizer_flat(self):
        """
        Test Ionizer on molecules without 3D coordinates.
        """
        self.mol.RemoveAllConformers()
        mol = Chem.RemoveHs(self.mol)
        assert self.mol.GetNumConformers() == 0
        self.ionized_mol.RemoveAllConformers()
        ref_mol = Chem.RemoveHs(self.ionized_mol)

        # make sure it calls the right method
        assert self.ionizer(mol).ToBinary() == self.ionizer._ionize_2d(
            mol).ToBinary()

        # compare molecule SMILES
        ionized_mol = self.ionizer(mol)
        assert Chem.MolToSmiles(
            ionized_mol, isomericSmiles=True) == Chem.MolToSmiles(
                ref_mol, isomericSmiles=True)

    def test_ionizer_conformers(self):
        """
        Make sure ionization preserves heavy atom coordinates.
        """
        assert self.mol.GetNumConformers() > 0

        # make sure it calls the right method
        assert self.ionizer(self.mol).ToBinary() == self.ionizer._ionize_3d(
            self.mol).ToBinary()

        # compare atomic coordinates
        ionized_mol = self.ionizer(self.mol)
        assert ionized_mol.GetNumConformers() > 0
        assert self.ionized_mol.GetNumConformers() > 0
        for a, b in zip(ionized_mol.GetConformers(),
                        self.ionized_mol.GetConformers()):
            for atom in ionized_mol.GetAtoms():
                if atom.GetAtomicNum() == 1:
                    continue  # hydrogens tend to move
                idx = atom.GetIdx()
                assert np.allclose(
                    list(a.GetAtomPosition(idx)), list(b.GetAtomPosition(idx)),
                    atol=0.00009)  # obabel rounds to four digits


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
