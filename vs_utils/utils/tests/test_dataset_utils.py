"""
Tests for dataset_utils.
"""
import gzip
import shutil
import tempfile
import unittest

from rdkit import Chem

from vs_utils.utils.dataset_utils import MoleculeDatabase


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

    def check_database(self, mols=None, filename=None):
        """
        Check database contents.

        Parameters
        ----------
        mols : list, optional
            Molecules that should be in the database. Defaults to self.mols.
        filename : str, optional
            Existing database filename.
        """
        if mols is None:
            mols = self.mols
        if filename is not None:
            database = MoleculeDatabase()
            database.load(filename)
        else:
            database = self.database

        # check for appropriate length
        assert len(database) == len(mols)

        # check that SMILES are what we expect
        for mol in mols:
            assert database.engine.get_smiles(mol) in database

    def test_load(self):
        """
        Test MoleculeDatabase.load.
        """
        _, filename = tempfile.mkstemp(dir=self.temp_dir)
        with open(filename, 'wb') as f:
            f.write('{}\n'.format(
                self.database.engine.get_smiles(self.mols[0])))
        self.check_database([self.mols[0]], filename)

    def test_load_gz(self):
        """
        Test MoleculeDatabase.load with gzipped input.
        """
        _, filename = tempfile.mkstemp(dir=self.temp_dir, suffix='.gz')
        with gzip.open(filename, 'wb') as f:
            f.write('{}\n'.format(
                self.database.engine.get_smiles(self.mols[0])))
        self.check_database([self.mols[0]], filename)

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
        self.check_database(filename=filename)
        with open(filename) as f:
            assert len(f.readlines()) == len(self.mols)

    def test_save_gz(self):
        """
        Test MoleculeDatabase.save with gzipped output.
        """
        _, filename = tempfile.mkstemp(dir=self.temp_dir, suffix='.gz')
        for mol in self.mols:
            self.database.add_mol(mol)
        self.database.save(filename)
        self.check_database(filename=filename)
        with gzip.open(filename) as f:
            assert len(f.readlines()) == len(self.mols)

    def test_add_mol(self):
        """
        Test MoleculeDatabase.add_mol.
        """
        for mol in self.mols:
            self.database.add_mol(mol)
        self.check_database()

    def test_add_mol_duplicate(self):
        """
        Test MoleculeDatabase.add_mol with a duplicate molecule.
        """
        for mol in self.mols:  # add once
            self.database.add_mol(mol)
        for mol in self.mols:  # add twice
            self.database.add_mol(mol)
        self.check_database()
