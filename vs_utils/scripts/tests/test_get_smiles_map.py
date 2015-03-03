"""
Tests for get_smiles_map.py.
"""
import shutil
import tempfile
import unittest

from rdkit import Chem

from vs_utils.scripts.get_smiles_map import main, parse_args
from vs_utils.utils import read_pickle


class TestGetSmilesMap(unittest.TestCase):
    """
    Tests for get_smiles_map.py.
    """
    def setUp(self):
        """
        Set up tests.
        """
        self.temp_dir = tempfile.mkdtemp()
        self.smiles = [
            'CC(=O)OC1=CC=CC=C1C(=O)O', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
            'CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F']
        self.cids = [2244, 3672, 2662]
        _, self.input_filename = tempfile.mkstemp(dir=self.temp_dir,
                                                  suffix='.smi')
        _, self.output_filename = tempfile.mkstemp(dir=self.temp_dir,
                                                   suffix='.pkl')

        # write SMILES to file
        with open(self.input_filename, 'wb') as f:
            for smile, cid in zip(self.smiles, self.cids):
                f.write('{}\t{}\n'.format(smile, cid))

    def tearDown(self):
        """
        Clean up tests.
        """
        shutil.rmtree(self.temp_dir)

    def test_main(self):
        """
        Test main.
        """
        args = parse_args(['-i', self.input_filename, '-o',
                           self.output_filename, '-p', 'CID'])
        main(args.input, args.output, args.prefix)
        data = read_pickle(self.output_filename)
        assert len(data) == len(self.smiles)
        for smile, cid in zip(self.smiles, self.cids):
            assert data['CID{}'.format(cid)] == Chem.MolToSmiles(
                Chem.MolFromSmiles(smile), isomericSmiles=True)

    def test_failure_on_bare_id(self):
        """
        Test failure on bare IDs.
        """
        args = parse_args(['-i', self.input_filename, '-o',
                           self.output_filename])
        try:
            main(args.input, args.output, args.prefix)
            raise AssertionError
        except TypeError:
            pass

    def test_update(self):
        """
        Test update existing map.
        """
        args = parse_args(['-i', self.input_filename, '-o',
                           self.output_filename, '-p', 'CID'])
        main(args.input, args.output, args.prefix, args.update)

        # add another molecule
        self.smiles.append('CC(=O)NC1=CC=C(C=C1)O')
        self.cids.append(1983)
        with open(self.input_filename, 'wb') as f:
            for smile, cid in zip(self.smiles, self.cids):
                f.write('{}\t{}\n'.format(smile, cid))

        # update existing map
        main(args.input, args.output, args.prefix, True)
        data = read_pickle(self.output_filename)
        assert len(data) == len(self.smiles)
        for smile, cid in zip(self.smiles, self.cids):
            assert data['CID{}'.format(cid)] == Chem.MolToSmiles(
                Chem.MolFromSmiles(smile), isomericSmiles=True)
