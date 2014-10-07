"""
Tests for target_utils.
"""
import numpy as np
import os
import shutil
import tempfile
import unittest

from .. import read_pickle, write_pickle
from ..target_utils import PcbaParser


class TestPcbaParser(unittest.TestCase):
    """
    Tests for PcbaParser.
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
        write_pickle(self.map, self.map_filename)

        # use a subset of AID588342
        # note that CID 654924 is duplicated
        this_dir = os.path.split(os.path.realpath(__file__))[0]
        self.data_filename = os.path.join(this_dir, 'data/test_assay_data.csv')

        # set up parser
        self.engine = PcbaParser(self.data_filename, self.map_filename)

    def tearDown(self):
        """
        Clean up tests.
        """
        shutil.rmtree(self.temp_dir)

    def test_read_pcba_data(self):
        """
        Test PcbaParser.read_pcba_data.
        """
        data = self.engine.read_pcba_data(self.data_filename)
        assert data.shape[0] == 4  # CID 654924 is duplicated
        assert np.array_equal(data.PUBCHEM_CID,
                              [645443, 645449, 654924, 2997889])

    def test_map_cids_to_smiles(self):
        """
        Test PcbaParser.map_cids_to_smiles.
        """
        data = self.engine.read_pcba_data(self.data_filename)
        id_map = read_pickle(self.map_filename)
        smiles, indices = self.engine.map_cids_to_smiles(
            data.PUBCHEM_CID, id_map)
        assert len(smiles) == len(indices) == 2
        assert smiles[0] == self.map['CID645443']
        assert smiles[1] == self.map['CID2997889']
        assert np.array_equal(indices, [0, 3])

    def test_get_targets_classification(self):
        """
        Test PcbaParser.get_targets with classification.
        """
        smiles, targets = self.engine.get_targets()
        assert len(smiles) == len(targets) == 2
        idx = np.where(smiles == self.map['CID2997889'])[0][0]
        assert targets[idx]  # marked Active
        idx = np.where(smiles == self.map['CID645443'])[0][0]
        assert not targets[idx]  # marked Inactive

    def test_get_targets_regression(self):
        """
        Test PcbaParser.get_targets with regression.
        """
        columns = [7, 8, 12, 14, 15, 20, 22, 23, 24, 25, 26]
        self.engine = PcbaParser(self.data_filename, self.map_filename,
                                 cols=columns)
        smiles, targets = self.engine.get_targets()
        idx = np.where(smiles == self.map['CID2997889'])[0][0]
        assert not np.any(np.isnan(targets[idx]))
        idx = np.where(smiles == self.map['CID645443'])[0][0]
        assert np.any(np.isnan(targets[idx]))  # will have NaNs
