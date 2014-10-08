"""
Tests for molecule_database.py.
"""
import shutil
import tempfile
import unittest

from rdkit import Chem

from pande_gas.scripts.molecule_database import main, parse_args


class TestMoleculeDatabase(unittest.TestCase):
    """
    Tests for molecule_database.py.
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

    def tearDown(self):
        """
        Clean up tests.
        """
        shutil.rmtree(self.temp_dir)

    def check_output(self, input_args):
        """
        Run main and return the resulting database.

        Parameters
        ----------
        args : list
            Command-line arguments.
        """
        args = parse_args(input_args)
        

    def test_main(self):
        """
        Test main.
        """

