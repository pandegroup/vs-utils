"""
Utilities for parsing target files from different sources (e.g. PubChem
BioAssay).
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

from collections import OrderedDict
import gzip
import numpy as np
import pandas as pd
import warnings

from rdkit_utils import serial

from pande_gas.utils import read_pickle, SmilesGenerator


class AssayDataParser(object):
    """
    Parse assay data files.

    Parameters
    ----------
    data_filename : str
        Data filename.
    map_filename : str
        Compound ID->SMILES map filename.
    primary_key : str
        Name of column containing compound IDs.
    id_prefix : str, optional
        Prefix to prepend to compound IDs for mapping compound IDs to SMILES.
    activity_key : str, optional
        Name of column containing compound activity assignments. Must be
        provided if column_indices is None. If both activity_key and
        column_indices are set, column_indices will be used.
    activity_value : str, optional
        Value of positive class in activity_key. For example, 'Active' when
        parsing PubChem BioAssay data.
    column_indices : list, optional
        Data column indices to include. Must be provided if activity_key is
        None. If both activity_key and column_indices are set, column_indices
        will be used.
    delimiter : str, optional (default '\t')
        Delimiter to use when parsing data file.
    """
    def __init__(self, data_filename, map_filename, primary_key,
                 id_prefix=None, activity_key=None, activity_value=None,
                 inactivity_value=None, column_indices=None, delimiter='\t'):
        self.data_filename = data_filename
        self.map_filename = map_filename
        self.primary_key = primary_key
        self.id_prefix = id_prefix
        if activity_key is None and column_indices is None:
            raise ValueError(
                'One of activity_key or column_indices must be set.')
        if activity_key is not None and activity_value is None:
            raise ValueError(
                'You must set activity_value when using activity_key.')
        self.activity_key = activity_key
        self.activity_value = activity_value
        self.inactivity_value = inactivity_value
        if column_indices is not None:
            column_indices = np.asarray(column_indices, dtype=int)
        self.column_indices = column_indices
        self.delimiter = delimiter

    def get_targets(self):
        """
        Parse data file and return targets and corresponding SMILES.

        Procedure
        ---------
        1. Read data and get unique rows by compound ID.
        2. Map compound IDs to SMILES.
        3. Extract targets from data.
        """
        data = self.read_data()
        id_map = read_pickle(self.map_filename)

        # get compound SMILES from map
        # indices are for data rows successfully mapped to SMILES
        smiles, indices = self.map_ids_to_smiles(data[self.primary_key],
                                                 id_map)

        # get targets
        if self.column_indices is not None:
            targets = np.zeros((data.shape[0], len(self.column_indices)),
                               dtype=float)
            for i, idx in enumerate(self.column_indices):
                targets[:, i] = data[data.columns[idx]]
        else:
            # targets other than self.activity_value and self.inactivity_value
            # are assigned -1
            targets = -1 * np.ones_like(data[self.activity_key], dtype=int)
            pos = (data[self.activity_key] == self.activity_value).values
            neg = (data[self.activity_key] == self.inactivity_value).values
            assert pos.sum() and neg.sum(), (
                'Dataset must have both active and inactive results.')
            targets[pos] = 1
            targets[neg] = 0
        targets = targets[indices]  # reduce targets to matched structures
        return smiles, targets

    def read_data(self, **kwargs):
        """
        Read assay data file.

        Parameters
        ----------
        kwargs : dict, optional
            Keyword arguments for pd.read_table.
        """
        if self.data_filename.endswith('.gz'):
            with gzip.open(self.data_filename) as f:
                df = pd.read_table(f, sep=self.delimiter, **kwargs)
        else:
            df = pd.read_table(self.data_filename, sep=self.delimiter,
                               **kwargs)
        df = df.drop_duplicates(self.primary_key)  # remove duplicate IDs
        return df

    def map_ids_to_smiles(self, ids, id_map):
        """
        Look up SMILES for compound IDs in a compound ID->SMILES map.

        Parameters
        ----------
        ids : array_like
            List of compound IDs.
        id_map : dict
            Compound ID->SMILES map.
        """
        smiles = []
        indices = []
        for i, this_id in enumerate(ids):
            if np.isnan(this_id):
                continue
            try:
                this_id = int(this_id)  # CIDs are often read in as floats
            except ValueError:
                pass
            if self.id_prefix is not None:
                # no bare IDs allowed in maps
                this_id = '{}{}'.format(self.id_prefix, this_id)
            if this_id in id_map:
                smiles.append(id_map[this_id])
                indices.append(i)
        return np.asarray(smiles), np.asarray(indices)

    def get_column_names(self):
        """
        Get names of selected data columns.
        """
        if self.column_indices is None:
            return
        names = []
        for i in self.column_indices:
            names.append(self.read_data().columns[i])
        return names


class PcbaParser(AssayDataParser):
    """
    Parse PubChem BioAssay (PCBA) target files.

    Parameters
    ----------
    data_filename : str
        Data filename.
    map_filename : str
        Compound ID->SMILES map filename.
    primary_key : str, optional (default 'PUBCHEM_CID')
        Name of column containing compound IDs.
    id_prefix : str, optional (default 'CID')
        Prefix to prepend to compound IDs for mapping compound IDs to SMILES.
    activity_key : str, optional (default 'PUBCHEM_ACTIVITY_OUTCOME')
        Name of column containing compound activity assignments. Must be
        provided if column_indices is None. If both activity_key and
        column_indices are set, column_indices will be used.
    activity_value : str, optional (default 'Active')
        Value of positive class in activity_key. For example, 'Active' when
        parsing PubChem BioAssay data.
    column_indices : list, optional
        Data column indices to include. Must be provided if activity_key is
        None. If both activity_key and column_indices are set, column_indices
        will be used.
    delimiter : str, optional (default ',')
        Delimiter to use when parsing data file.
    """
    def __init__(self, data_filename, map_filename, primary_key='PUBCHEM_CID',
                 id_prefix='CID', activity_key='PUBCHEM_ACTIVITY_OUTCOME',
                 activity_value='Active', inactivity_value='Inactive',
                 column_indices=None, delimiter=','):
        super(PcbaParser, self).__init__(
            data_filename, map_filename, primary_key, id_prefix, activity_key,
            activity_value, inactivity_value, column_indices, delimiter)


class Nci60Parser(AssayDataParser):
    """
    Parse NCI60 target file.

    Parameters
    ----------
    data_filename : str
        Data filename.
    map_filename : str
        Compound ID->SMILES map filename.
    primary_key : str, optional (default 'NSC')
        Name of column containing compound IDs.
    id_prefix : str, optional (default 'NSC')
        Prefix to prepend to compound IDs for mapping compound IDs to SMILES.
    activity_key : str, optional
        Name of column containing compound activity assignments. Must be
        provided if column_indices is None. If both activity_key and
        column_indices are set, column_indices will be used.
    activity_value : str, optional
        Value of positive class in activity_key. For example, 'Active' when
        parsing PubChem BioAssay data.
    column_indices : list, optional (default range(4, 64))
        Data column indices to include. Must be provided if activity_key is
        None. If both activity_key and column_indices are set, column_indices
        will be used.
    delimiter : str, optional (default '\t')
        Delimiter to use when parsing data file.
    """
    def __init__(self, data_filename, map_filename, primary_key='NSC',
                 id_prefix='NSC', activity_key=None, activity_value=None,
                 inactivity_value=None, column_indices=range(4, 64),
                 delimiter='\t'):
        super(Nci60Parser, self).__init__(
            data_filename, map_filename, primary_key, id_prefix, activity_key,
            activity_value, inactivity_value, column_indices, delimiter)

    def read_data(self, **kwargs):
        """
        Read assay data file.

        Parameters
        ----------
        kwargs : dict, optional
            Keyword arguments for pd.read_table.
        """
        # treat '-' and 'na' values as NaNs
        return super(Nci60Parser, self).read_data(na_values=['-', 'na'])

    def split_targets(self):
        """
        Split targets among different assays.
        """
        df = self.read_data()
        names = df.columns[self.column_indices]
        smiles, targets = self.get_targets()
        split_targets = OrderedDict()
        for i, name in enumerate(names):
            keep = ~np.isnan(targets[:, i])
            if not np.count_nonzero(keep):
                warnings.warn(
                    'Assay "{}" has no matching records.'.format(name))
                continue
            split_targets[name] = {'smiles': smiles[keep],
                                   'targets': targets[keep]}
        return split_targets


class Tox21Parser(object):
    """
    Parse Tox21 data files.

    Parameters
    ----------
    filename : str
        Data filename.
    merge_strategy : str, optional (default 'max')
        Strategy to use when merging targets for duplicated molecules. Choose
        from 'max' (active if active in any assay), 'min' (inactive if inactive
        in any assay), 'majority_pos' (majority vote with ties assigned
        active), or 'majority_neg' (majority vote with ties assigned inactive).
    """
    dataset_names = ['NR-AR', 'NR-AhR', 'NR-AR-LBD', 'NR-ER', 'NR-ER-LBD',
                     'NR-Aromatase', 'NR-PPAR-gamma', 'SR-ARE', 'SR-ATAD5',
                     'SR-HSE', 'SR-MMP', 'SR-p53']

    def __init__(self, filename, merge_strategy='max'):
        self.filename = filename
        assert merge_strategy in ['max', 'min', 'majority_pos', 'majority_neg']
        self.merge_strategy = merge_strategy

    def read_data(self):
        """
        Read labeled molecules.
        """
        with serial.MolReader().open(self.filename) as reader:
            mols = list(reader)
        return mols

    def read_targets(self):
        """
        Get labels for molecules from SD data fields matching dataset names.

        Returns
        -------
        data : dict
            Nested dictionary containing SMILES and targets for compounds in
            each dataset. Keyed by data->dataset->SMILES->target, where target
            is a list.
        """
        engine = SmilesGenerator()
        data = {dataset: {} for dataset in self.dataset_names}
        skipped = []
        for mol in self.read_data():
            smiles = engine.get_smiles(mol)
            for prop in list(mol.GetPropNames()):
                if prop in data:
                    score = int(mol.GetProp(prop))
                    if smiles not in data[prop]:
                        data[prop][smiles] = []
                    data[prop][smiles].append(score)
                else:  # skip irrelevant SD fields
                    if prop not in skipped:
                        skipped.append(prop)
                    continue
        print 'Skipped properties:\n{}'.format('\n'.join(skipped))
        return data

    def merge_targets(self, data):
        """
        Merge labels for duplicate molecules according to a specified merge
        stratecy ('max', 'min', 'majority_pos', 'majority_neg').

        Parameters
        ----------
        data : dict
            Nested dictionary containing SMILES and targets for compounds in
            each dataset. Keyed by data->dataset->SMILES->target, where target
            is a list.

        Returns
        -------
        data : dict
            Nested dictionary containing SMILES and targets for compounds in
            each dataset. Keyed by data->dataset->SMILES->target, where target
            is an integer.
        """
        for dataset in self.dataset_names:
            for smiles, targets in data[dataset].items():
                targets = np.asarray(targets, dtype=int)
                if self.merge_strategy == 'max':
                    data[dataset][smiles] = max(targets)
                elif self.merge_strategy == 'min':
                    data[dataset][smiles] = min(targets)
                # 0.5 rounds down
                elif self.merge_strategy == 'majority_neg':
                    data[dataset][smiles] = int(np.round(np.mean(targets)))
                # 0.5 rounds up
                elif self.merge_strategy == 'majority_pos':
                    data[dataset][smiles] = (int(np.round(
                        np.mean(targets) + 1)) - 1)
        return data

    def get_targets(self):
        """
        Get SMILES and targets for each Tox21 dataset.
        """
        split_targets = {}
        data = self.merge_targets(self.read_targets())
        for dataset in data:
            if not len(data[dataset]):
                warnings.warn('Dataset "{}" is empty'.format(dataset))
                continue
            smiles, targets = [], []
            for this_smiles, target in data[dataset].items():
                smiles.append(this_smiles)
                targets.append(target)
            split_targets[dataset] = {'smiles': np.asarray(smiles),
                                      'targets': np.asarray(
                                          targets, dtype=int)}
        return split_targets


class Counterscreen(object):
    """
    Reassign labels based on counterscreen results.

    Parameters
    ----------
    primary : AssayDataParser
        Primary assay.
    counter : list
        Counterscreens.
    """
    def __init__(self, primary, counter):
        self.primary = primary
        self.counter = counter

    def get_targets(self):
        """
        Reassign labels based on counterscreen results. Compounds marked
        'Active' in primary and counter will be reassigned to label -2.
        """
        smiles, targets = self.primary.get_targets()
        for counter in self.counter:
            counter_smiles, counter_targets = counter.get_targets()
            counter_actives = counter_smiles[counter_targets == 1]
            targets[np.in1d(smiles, counter_actives)] = -2
        return smiles, targets