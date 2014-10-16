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

from pande_gas.utils import read_pickle


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
                 column_indices=None, delimiter='\t'):
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
            targets = np.asarray(
                data[self.activity_key] == self.activity_value)
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
                 activity_value='Active', column_indices=None, delimiter=','):
        super(PcbaParser, self).__init__(
            data_filename, map_filename, primary_key, id_prefix, activity_key,
            activity_value, column_indices, delimiter)


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
                 column_indices=range(4, 64), delimiter='\t'):
        super(Nci60Parser, self).__init__(
            data_filename, map_filename, primary_key, id_prefix, activity_key,
            activity_value, column_indices, delimiter)

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
