"""
Tests for target_utils.
"""
import numpy as np
import os
import shutil
import tempfile
import unittest

from rdkit import Chem

from rdkit_utils import serial

from .. import read_pickle, write_pickle
from ..target_utils import (AssayDataParser, Counterscreen, Nci60Parser,
                            PcbaParser, Tox21Parser)


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
            'CID654924': 'Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1',
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
                                      activity_value='Active',
                                      inactivity_value='Inactive',
                                      id_prefix='CID')

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
        assert len(smiles) == len(indices) == 3
        assert smiles[0] == self.map['CID645443']
        assert smiles[1] == self.map['CID654924']
        assert smiles[2] == self.map['CID2997889']
        assert np.array_equal(indices, [0, 2, 3]), indices

    def test_inconclusive(self):
        """
        Test labeling of 'Inconclusive' results.
        """
        smiles, targets = self.engine.get_targets()
        assert np.array_equal(targets, [0, -1, 1]), targets


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


class TestTox21Parser(unittest.TestCase):
    """
    Test Tox21Parser.
    """
    def setUp(self):
        """
        Set up tests.
        """
        self.temp_dir = tempfile.mkdtemp()
        self.mols = []
        aspirin = Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(=O)O')
        aspirin.SetProp('_Name', 'aspirin1')
        aspirin.SetProp('NR-Aromatase', '0')
        aspirin.SetProp('SR-ATAD5', '1')
        aspirin.SetProp('NR-ER', '0')
        aspirin.SetProp('SR-p53', '0')
        aspirin.SetProp('NR-AhR', '1')
        aspirin.SetProp('SR-ARE', '1')
        aspirin.SetProp('Irrelevant', 'blah')
        self.mols.append(aspirin)
        ibuprofen = Chem.MolFromSmiles('CC(C)CC1=CC=C(C=C1)C(C)C(=O)O')
        ibuprofen.SetProp('_Name', 'ibuprofen')
        ibuprofen.SetProp('SR-ATAD5', '1')
        ibuprofen.SetProp('NR-ER', '1')
        ibuprofen.SetProp('SR-p53', '1')
        ibuprofen.SetProp('NR-AhR', '1')
        ibuprofen.SetProp('SR-ARE', '0')
        ibuprofen.SetProp('Irrelevant', 'blah2')
        self.mols.append(ibuprofen)
        aspirin2 = Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(=O)O')
        aspirin2.SetProp('_Name', 'aspirin2')
        aspirin2.SetProp('SR-ATAD5', '0')
        aspirin2.SetProp('NR-ER', '0')
        aspirin2.SetProp('SR-p53', '0')
        aspirin2.SetProp('NR-AhR', '1')
        aspirin2.SetProp('SR-ARE', '1')
        aspirin2.SetProp('Irrelevant', 'blah3')
        self.mols.append(aspirin2)
        aspirin3 = Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(=O)O')
        aspirin3.SetProp('_Name', 'aspirin3')
        aspirin3.SetProp('SR-ATAD5', '1')
        aspirin3.SetProp('NR-ER', '0')
        aspirin3.SetProp('SR-p53', '0')
        aspirin3.SetProp('NR-AhR', '1')
        aspirin3.SetProp('SR-ARE', '1')
        aspirin3.SetProp('Irrelevant', 'blah4')
        self.mols.append(aspirin3)
        aspirin4 = Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(=O)O')
        aspirin4.SetProp('_Name', 'aspirin4')
        aspirin4.SetProp('SR-ATAD5', '0')
        aspirin4.SetProp('NR-ER', '0')
        aspirin4.SetProp('SR-p53', '1')
        aspirin4.SetProp('NR-AhR', '0')
        aspirin4.SetProp('SR-ARE', '1')
        aspirin4.SetProp('Irrelevant', 'blah5')
        self.mols.append(aspirin4)
        self.smiles = [Chem.MolToSmiles(mol, isomericSmiles=True)
                       for mol in self.mols]

        # write input file
        _, self.filename = tempfile.mkstemp(dir=self.temp_dir, suffix='.sdf')
        with serial.MolWriter().open(self.filename) as writer:
            writer.write(self.mols)

        # set up parser
        self.engine = Tox21Parser(self.filename)

    def tearDown(self):
        """
        Clean up tests.
        """
        shutil.rmtree(self.temp_dir)

    def test_read_data(self):
        """
        Test Tox21Parser.read_data.
        """
        assert len(self.engine.read_data()) == 5

    def test_read_targets(self):
        """
        Test Tox21Parser.read_targets.
        """
        data = self.engine.read_targets()
        count = 0
        for dataset in data:
            assert dataset in self.engine.dataset_names
            if len(data[dataset]):
                count += 1
            else:
                continue
        assert count == 6

        # check individual datasets
        for dataset in data:
            if not len(data[dataset]):
                continue
            if dataset == 'NR-Aromatase':
                assert set(data[dataset].keys()) == {self.smiles[0]}
            else:
                assert set(data[dataset].keys()) == set(self.smiles)

        assert np.array_equal(data['NR-Aromatase'].values(), [[0]])

        assert [1, 0, 1, 0] in data['SR-ATAD5'].values()
        assert [1] in data['SR-ATAD5'].values()
        assert [0] not in data['SR-ATAD5'].values()

        assert [0, 0, 0, 0] in data['NR-ER'].values()
        assert [1] in data['NR-ER'].values()
        assert [0] not in data['NR-ER'].values()

        assert [0, 0, 0, 1] in data['SR-p53'].values()
        assert [1] in data['SR-p53'].values()
        assert [0] not in data['SR-p53'].values()

        assert [1, 1, 1, 0] in data['NR-AhR'].values()
        assert [1] in data['NR-AhR'].values()
        assert [0] not in data['NR-AhR'].values()

        assert [1, 1, 1, 1] in data['SR-ARE'].values()
        assert [0] in data['SR-ARE'].values()
        assert [1] not in data['SR-ARE'].values()

    def test_merge_targets_max(self):
        """
        Test Tox21Parser.merge_targets with 'max' merge_strategy.
        """
        self.engine.merge_strategy = 'max'
        data = self.engine.merge_targets(self.engine.read_targets())

        assert np.array_equal(data['NR-Aromatase'].values(), [0])

        assert data['SR-ATAD5'][self.smiles[0]] == 1
        assert data['SR-ATAD5'][self.smiles[1]] == 1

        assert data['NR-ER'][self.smiles[0]] == 0
        assert data['NR-ER'][self.smiles[1]] == 1

        assert data['SR-p53'][self.smiles[0]] == 1
        assert data['SR-p53'][self.smiles[1]] == 1

        assert data['NR-AhR'][self.smiles[0]] == 1
        assert data['NR-AhR'][self.smiles[1]] == 1

        assert data['SR-ARE'][self.smiles[0]] == 1
        assert data['SR-ARE'][self.smiles[1]] == 0

    def test_merge_targets_min(self):
        """
        Test Tox21Parser.merge_targets with 'min' merge_strategy.
        """
        self.engine.merge_strategy = 'min'
        data = self.engine.merge_targets(self.engine.read_targets())

        assert np.array_equal(data['NR-Aromatase'].values(), [0])

        assert data['SR-ATAD5'][self.smiles[0]] == 0
        assert data['SR-ATAD5'][self.smiles[1]] == 1

        assert data['NR-ER'][self.smiles[0]] == 0
        assert data['NR-ER'][self.smiles[1]] == 1

        assert data['SR-p53'][self.smiles[0]] == 0
        assert data['SR-p53'][self.smiles[1]] == 1

        assert data['NR-AhR'][self.smiles[0]] == 0
        assert data['NR-AhR'][self.smiles[1]] == 1

        assert data['SR-ARE'][self.smiles[0]] == 1
        assert data['SR-ARE'][self.smiles[1]] == 0

    def test_merge_targets_majority_pos(self):
        """
        Test Tox21Parser.merge_targets with 'majority_pos' merge_strategy.
        """
        self.engine.merge_strategy = 'majority_pos'
        data = self.engine.merge_targets(self.engine.read_targets())

        assert np.array_equal(data['NR-Aromatase'].values(), [0])

        assert data['SR-ATAD5'][self.smiles[0]] == 1
        assert data['SR-ATAD5'][self.smiles[1]] == 1

        assert data['NR-ER'][self.smiles[0]] == 0
        assert data['NR-ER'][self.smiles[1]] == 1

        assert data['SR-p53'][self.smiles[0]] == 0
        assert data['SR-p53'][self.smiles[1]] == 1

        assert data['NR-AhR'][self.smiles[0]] == 1
        assert data['NR-AhR'][self.smiles[1]] == 1

        assert data['SR-ARE'][self.smiles[0]] == 1
        assert data['SR-ARE'][self.smiles[1]] == 0

    def test_merge_targets_majority_neg(self):
        """
        Test Tox21Parser.merge_targets with 'majority_neg' merge_strategy.
        """
        self.engine.merge_strategy = 'majority_neg'
        data = self.engine.merge_targets(self.engine.read_targets())

        assert np.array_equal(data['NR-Aromatase'].values(), [0])

        assert data['SR-ATAD5'][self.smiles[0]] == 0
        assert data['SR-ATAD5'][self.smiles[1]] == 1

        assert data['NR-ER'][self.smiles[0]] == 0
        assert data['NR-ER'][self.smiles[1]] == 1

        assert data['SR-p53'][self.smiles[0]] == 0
        assert data['SR-p53'][self.smiles[1]] == 1

        assert data['NR-AhR'][self.smiles[0]] == 1
        assert data['NR-AhR'][self.smiles[1]] == 1

        assert data['SR-ARE'][self.smiles[0]] == 1
        assert data['SR-ARE'][self.smiles[1]] == 0

    def test_get_targets(self):
        """
        Test Tox21Parser.get_targets.
        """
        data = self.engine.get_targets()
        assert len(data) == 6
        for dataset in data:
            smiles = data[dataset]['smiles']
            targets = data[dataset]['targets']
            assert len(smiles) == len(targets) > 0


class TestCounterscreen(TestAssayDataParser):
    """
    Tests for Counterscreen.
    """
    def test_counterscreen(self):
        """
        Test reassignment of counterscreen actives.
        """
        counter = Counterscreen(self.engine, [self.engine])
        smiles, targets = counter.get_targets()
        assert np.array_equal(targets, [0, -1, -2]), targets
