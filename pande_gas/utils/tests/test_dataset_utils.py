"""
Tests for dataset_utils.
"""
import shutil
import tempfile
import unittest

from rdkit import Chem

from ..dataset_utils import MoleculeDatabase


class TestMoleculeDatabase(unittest.TestCase):
    """
    Tests for MoleculeDatabase.
    """
    def setUp(self):
        """
        Set up tests.
        """
        smiles = ['CC(=O)OC1=CC=CC=C1C(=O)O', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
                  'CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F']
        names = ['aspirin', 'ibuprofen', 'celecoxib']
        self.cids = [2244, 3672, 2662]
        self.mols = []
        for s, n in zip(smiles, names):
            mol = Chem.MolFromSmiles(s)
            mol.SetProp('_Name', n)
            self.mols.append(mol)
        self.temp_dir = tempfile.mkdtemp()
        self.database = MoleculeDatabase()

    def tearDown(self):
        """
        Clean up tests.
        """
        shutil.rmtree(self.temp_dir)

    def test_load(self):
        """
        Test MoleculeDatabase.load.
        """
        _, filename = tempfile.mkstemp(dir=self.temp_dir)
        with open(filename, 'wb') as f:
            f.write('{}\n'.format(
                self.database.engine.get_smiles(self.mols[0])))
        self.database.load(filename)
        assert len(self.database) == 1
        for mol in self.mols:
            self.database.add_mol(mol)
        assert len(self.database) == len(self.mols)  # try adding duplicates

    def test_load_bogus(self):
        """
        Test failure on loading a bogus dataset.
        """
        _, filename = tempfile.mkstemp(dir=self.temp_dir)
        with open(filename, 'wb') as f:
            f.write('bogus\n')
        try:
            self.database.load(filename)
            raise AssertionError
        except ValueError:
            pass

    def test_save(self):
        """
        Test MoleculeDatabase.save.
        """
        _, filename = tempfile.mkstemp(dir=self.temp_dir)
        for mol in self.mols:
            self.database.add_mol(mol)
        self.database.save(filename)
        self.database = MoleculeDatabase()
        self.database.load(filename)
        assert len(self.database) == len(self.mols)

    def test_add_mol(self):
        """
        Test MoleculeDatabase.add_mol.
        """
        for mol in self.mols:
            self.database.add_mol(mol)
        assert len(self.database) == len(self.mols)

    def test_add_mol_duplicate(self):
        """
        Test MoleculeDatabase.add_mol with a duplicate molecule.
        """
        for mol in self.mols:  # add once
            self.database.add_mol(mol)
        for mol in self.mols:  # add twice
            self.database.add_mol(mol)
        assert len(self.database) == len(self.mols)
