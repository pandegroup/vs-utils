"""
Tests for parse_tox21_datasets.py.
"""
import glob
import os

from rdkit import Chem

from vs_utils.scripts.datasets.parse_tox21_datasets import main, get_args
from vs_utils.utils import read_pickle
from vs_utils.utils.tests.test_target_utils import TestTox21Parser


class TestParseTox21Datasets(TestTox21Parser):
    """
    Tests for parse_tox21_datasets.py.
    """
    def test_main(self):
        """
        Test main.
        """
        args = get_args([self.filename, '-d', self.temp_dir])
        main(args.input, args.merge, args.dir)

        # check for the right number of files
        assert len(glob.glob(os.path.join(self.temp_dir, '*.pkl.gz'))) == 6

        # inspect files individually
        for filename in glob.glob(os.path.join(self.temp_dir, '*.pkl.gz')):
            data = read_pickle(filename)
            assert len(data['smiles']) == len(data['targets'])

            # try to read SMILES
            for this_smiles in data['smiles']:
                Chem.MolFromSmiles(this_smiles)

            # check type of targets
            assert data['targets'].dtype == int
