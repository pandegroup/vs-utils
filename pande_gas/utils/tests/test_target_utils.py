"""
Tests for target_utils.
"""
import numpy as np
import os
import shutil
import tempfile
import unittest

from .. import read_pickle, write_pickle
from ..target_utils import AssayDataParser, Nci60Parser, PcbaParser


class TestAssayDataParser(unittest.TestCase):
    """
    Tests for AssayDataParser.
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
        self.data_filename = os.path.join(this_dir, 'data/test_pcba_data.csv')

        # set up parser
        # settings match PcbaParser defaults
        self.engine = AssayDataParser(self.data_filename, self.map_filename,
                                      delimiter=',', primary_key='PUBCHEM_CID',
                                      activity_key='PUBCHEM_ACTIVITY_OUTCOME',
                                      activity_value='Active', id_prefix='CID')

    def tearDown(self):
        """
        Clean up tests.
        """
        shutil.rmtree(self.temp_dir)

    def test_read_data(self):
        """
        Test AssayDataParser.read_data.
        """
        data = self.engine.read_data()
        assert data.shape[0] == 4  # CID 654924 is duplicated
        assert np.array_equal(data.PUBCHEM_CID,
                              [645443, 645449, 654924, 2997889])

    def test_map_ids_to_smiles(self):
        """
        Test AssayDataParser.map_ids_to_smiles.
        """
        data = self.engine.read_data()
        id_map = read_pickle(self.map_filename)
        smiles, indices = self.engine.map_ids_to_smiles(
            data.PUBCHEM_CID, id_map)
        assert len(smiles) == len(indices) == 2
        assert smiles[0] == self.map['CID645443']
        assert smiles[1] == self.map['CID2997889']
        assert np.array_equal(indices, [0, 3])


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
        self.data_filename = os.path.join(this_dir, 'data/test_pcba_data.csv')

        # set up parser
        self.engine = PcbaParser(self.data_filename, self.map_filename)

    def tearDown(self):
        """
        Clean up tests.
        """
        shutil.rmtree(self.temp_dir)

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
                                 column_indices=columns)
        smiles, targets = self.engine.get_targets()
        idx = np.where(smiles == self.map['CID2997889'])[0][0]
        assert not np.any(np.isnan(targets[idx]))
        idx = np.where(smiles == self.map['CID645443'])[0][0]
        assert np.any(np.isnan(targets[idx]))  # will have NaNs


class TestNci60Parser(unittest.TestCase):
    """
    Tests for Nci60Parser.
    """
    def setUp(self):
        """
        Set up tests.
        """
        self.temp_dir = tempfile.mkdtemp()
        self.map = {'NSC1': 'CC1=CC(=O)C=CC1=O'}
        _, self.map_filename = tempfile.mkstemp(dir=self.temp_dir,
                                                suffix='.pkl')
        write_pickle(self.map, self.map_filename)

        this_dir = os.path.split(os.path.realpath(__file__))[0]
        self.data_filename = os.path.join(this_dir, 'data/test_nci60_data.txt')

        # set up parser
        self.engine = Nci60Parser(self.data_filename, self.map_filename)

    def tearDown(self):
        """
        Clean up tests.
        """
        shutil.rmtree(self.temp_dir)

    def test_read_data(self):
        """
        Test Nci60Parser.read_data.
        """
        df = self.engine.read_data()
        fixed_count = df.count().values.sum()  # count excluding NaNs
        # use PcbaParser to read data (w/o proper NaN handling)
        engine = PcbaParser(self.data_filename, self.map_filename,
                            delimiter='\t', primary_key='NSC', id_prefix='NSC')
        df = engine.read_data()
        broken_count = df.count().values.sum()
        assert fixed_count < broken_count

    def test_get_targets(self):
        """
        Test Nci60Parser.get_targets.
        """
        smiles, targets = self.engine.get_targets()
        assert len(smiles) == len(targets) == 1
        assert smiles[0] == self.map.values()[0]
        assert targets.shape == (1, 60)

    def test_split_targets(self):
        """
        Test Nci60Parser.split_targets.
        """
        split_targets = self.engine.split_targets()
        assert len(split_targets) == 59  # one is empty
        for name in self.engine.get_column_names():
            if name == 'ME:MDA_N':
                continue  # this assay is empty
            assert name in split_targets
