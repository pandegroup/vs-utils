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
        self.mol = Chem.MolFromSmiles(smiles)

        # ionize the oxygen manually
        self.ionized_mol = Chem.Mol(self.mol)  # create a copy
        found = False
        for atom in self.ionized_mol.GetAtoms():
            if atom.GetAtomicNum() == 8 and atom.GetNumImplicitHs() == 1:
                atom.SetFormalCharge(-1)
                atom.SetNoImplicit(True)
                found = True
        assert found

        self.ionizer = ob_utils.Ionizer()

    def test_ionizer_flat(self):
        """
        Test Ionizer on molecules without 3D coordinates.
        """
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

        # generate conformers
        engine = conformers.ConformerGenerator(max_conformers=1)
        mol = engine.generate_conformers(self.mol)
        assert mol.GetNumConformers() > 0
        ref_mol = engine.generate_conformers(self.ionized_mol)
        assert ref_mol.GetNumConformers() > 0
        ref_mol = Chem.RemoveHs(ref_mol)

        # make sure Ionizer calls the right method
        assert self.ionizer(mol).ToBinary() != self.ionizer._ionize_2d(
            mol).ToBinary()
        assert self.ionizer(mol).ToBinary() == self.ionizer._ionize_3d(
            mol).ToBinary()

        # compare molecule SMILES
        ionized_mol = self.ionizer(mol)
        ionized_mol = Chem.RemoveHs(ionized_mol)
        assert Chem.MolToSmiles(
            ionized_mol, isomericSmiles=True) == Chem.MolToSmiles(
                self.ionized_mol, isomericSmiles=True)

        # compare heavy atom coordinates
        assert ionized_mol.GetNumConformers() > 0
        assert ionized_mol.GetNumConformers() == ref_mol.GetNumConformers()
        for a, b in zip(ionized_mol.GetConformers(),
                        ref_mol.GetConformers()):
            for atom in ionized_mol.GetAtoms():
                if atom.GetAtomicNum() == 1:
                    continue  # hydrogens tend to move
                idx = atom.GetIdx()
                a_pos = list(a.GetAtomPosition(idx))
                b_pos = list(b.GetAtomPosition(idx))

                # obabel rounds to four digits
                assert np.allclose(a_pos, b_pos, atol=0.00009), (a_pos, b_pos)


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
