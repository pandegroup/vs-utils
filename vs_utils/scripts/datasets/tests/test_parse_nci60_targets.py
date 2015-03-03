"""
Tests for parse_nci60_targets.py.
"""
import cPickle
import glob
import gzip
import os
import shutil
import tempfile
import unittest

from vs_utils.scripts.datasets.parse_nci60_targets import main, parse_args


class TestParseNci60Targets(unittest.TestCase):
    """
    Tests for parse_pcba_targets.py.
    """
    def setUp(self):
        """
        Set up tests.
        """
        self.temp_dir = tempfile.mkdtemp()
        self.map = {'NSC1': 'CC1=CC(=O)C=CC1=O'}
        _, self.map_filename = tempfile.mkstemp(dir=self.temp_dir,
                                                suffix='.pkl')
        with open(self.map_filename, 'wb') as f:
            cPickle.dump(self.map, f, cPickle.HIGHEST_PROTOCOL)

        this_dir = os.path.split(os.path.realpath(__file__))[0]
        self.data_filename = os.path.join(
            this_dir, '../../../utils/tests/data/test_nci60_data.txt')

    def tearDown(self):
        """
        Clean up tests.
        """
        shutil.rmtree(self.temp_dir)

    def test_output(self):
        """
        Check output.
        """
        args = parse_args(['-i', self.data_filename, '-m', self.map_filename,
                           '-d', self.temp_dir])
        main(args.input, args.map, args.dir, args.prefix, args.suffix)
        count = 0
        for filename in glob.glob('{}/*.pkl.gz'.format(self.temp_dir)):
            count += 1
            with gzip.open(filename) as f:
                data = cPickle.load(f)
            assert len(data['smiles']) == len(data['targets'])
            assert len(data['targets']) == 1, len(data['targets'])
        assert count == 59, count  # one target is NaN
