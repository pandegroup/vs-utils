"""
Tests for parse_pcba_targets.py.
"""
import cPickle
import numpy as np
import os
import shutil
import tempfile
import unittest

from pande_gas.scripts.parse_pcba_targets import main, parse_args


class TestParsePcbaTargets(unittest.TestCase):
    """
    Tests for parse_pcba_targets.py.
    """
    def setUp(self):
        """
        Set up tests.
        """
        self.temp_dir = tempfile.mkdtemp()
        self.map = {
            'CID645443':
            'Cc1ccc(-n2c3c(cc(C(=O)Nc4cccc(C)n4)c2=O)C(=O)CCC3)cc1',
            'CID2997889': 'CC(C)(C)C(=O)Nc1ccc(-c2cn3ccsc3n2)cc1',
            'CID2244': 'CC(=O)Oc1ccccc1C(=O)O',
            'CID2662': 'Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1',
            'CID3672': 'CC(C)Cc1ccc(C(C)C(=O)O)cc1'}
        _, self.map_filename = tempfile.mkstemp(dir=self.temp_dir,
                                                suffix='.pkl')
        with open(self.map_filename, 'wb') as f:
            cPickle.dump(self.map, f, cPickle.HIGHEST_PROTOCOL)

        # use a subset of AID588342
        this_dir = os.path.split(os.path.realpath(__file__))[0]
        self.data_filename = os.path.join(this_dir, 'data/test_assay_data.csv')
        _, self.output_filename = tempfile.mkstemp(dir=self.temp_dir,
                                                   suffix='.pkl')

    def tearDown(self):
        """
        Clean up tests.
        """
        shutil.rmtree(self.temp_dir)

    def run_script(self, input_args):
        """
        Run script with the given arguments.

        Parameters
        ----------
        input_args : list
            Arguments.
        """
        args = parse_args(input_args)
        main(args.input, args.map, args.output, args.cols)
        with open(self.output_filename) as f:
            data = cPickle.load(f)
        assert len(data['smiles']) == len(data['targets']) == 2
        return data

    def test_classification_main(self):
        """
        Test classification.
        """
        args = ['-i', self.data_filename, '-m', self.map_filename, '-o',
                self.output_filename]
        data = self.run_script(args)
        idx = np.where(data['smiles'] == self.map['CID2997889'])[0][0]
        assert data['targets'][idx]  # marked Active
        idx = np.where(data['smiles'] == self.map['CID645443'])[0][0]
        assert not data['targets'][idx]  # marked Inactive

    def test_regression_main(self):
        """
        Test regression.
        """
        columns = ['7', '8', '12', '14', '15', '20', '22', '23', '24', '25',
                   '26']
        args = ['-i', self.data_filename, '-m', self.map_filename, '-o',
                self.output_filename, '-c'] + columns
        data = self.run_script(args)
        idx = np.where(data['smiles'] == self.map['CID2997889'])[0][0]
        assert not np.any(np.isnan(data['targets'][idx]))
        idx = np.where(data['smiles'] == self.map['CID645443'])[0][0]
        assert np.any(np.isnan(data['targets'][idx]))  # will have NaNs
