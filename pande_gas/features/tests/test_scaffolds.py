"""
Tests for scaffolds.py.
"""
import unittest

from rdkit import Chem

from ..scaffolds import Scaffold


class TestScaffold(unittest.TestCase):
    """
    Test Scaffold.
    """
    def setUp(self):
        """
        Set up tests.
        """
        smiles = ['CC(=O)OC1=CC=CC=C1C(=O)O',
                  'CN1C=C(C2=CC=CC=C21)C(=O)[C@@H]3CCC4=C(C3)NC=N4']
        names = ['aspirin', 'ramosetron']
        self.mols = []
        for this_smiles, name in zip(smiles, names):
            mol = Chem.MolFromSmiles(this_smiles)
            mol.SetProp('_Name', name)
            self.mols.append(mol)
        self.engine = Scaffold()

    def test_scaffolds(self):
        """
        Test scaffold generation.
        """
        scaffolds = self.engine(self.mols)
        assert len(scaffolds) == 2
        scaffold_mols = [Chem.MolFromSmiles(scaffold)
                         for scaffold in scaffolds]
        for mol, ref_mol in zip(scaffold_mols, self.mols):
            assert mol.GetNumAtoms() < ref_mol.GetNumAtoms()
        assert scaffold_mols[0].GetNumAtoms() == 6
        assert scaffold_mols[1].GetNumAtoms() == 20

    def test_chiral_scaffolds(self):
        """
        Test chiral scaffold generation.
        """
        achiral_scaffold = self.engine([self.mols[1]])[0]
        self.engine = Scaffold(include_chirality=True)
        chiral_scaffold = self.engine([self.mols[1]])[0]
        assert '@' not in achiral_scaffold
        assert '@' in chiral_scaffold
        assert (Chem.MolFromSmiles(achiral_scaffold).GetNumAtoms() ==
                Chem.MolFromSmiles(chiral_scaffold).GetNumAtoms())
