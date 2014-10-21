"""
Tests for parse_tox21_datasets.py.
"""
import glob

from rdkit import Chem

from pande_gas.scripts.datasets.parse_tox21_datasets import main, get_args
from pande_gas.utils import read_pickle
from pande_gas.utils.tests.test_target_utils import TestTox21Parser


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

        count = 0
        for filename in glob.glob('*.pkl.gz'):
            count += 1
            data = read_pickle(filename)
            assert len(data['smiles']) == len(data['targets'])

            # try to read SMILES
            for this_smiles in data['smiles']:
                Chem.MolFromSmiles(this_smiles)

            # check type of targets
            assert data['targets'].dtype == int
        assert count == 6
