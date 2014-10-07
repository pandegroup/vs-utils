"""
Utilities for parsing target files from different sources (e.g. PubChem
BioAssay).
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import gzip
import numpy as np
import pandas as pd

from pande_gas.utils import read_pickle


class PcbaParser(object):
    """
    Parse PubChem BioAssay (PCBA) target files.

    Parameters
    ----------
    data_filename : str
        PCBA data filename.
    map_filename : str
        ID->SMILES map filename.
    cols : list, optional
        Data column indices to include. If None, compounds are classified by
        assigned activity designations.
    """
    def __init__(self, data_filename, map_filename, cols=None):
        self.data_filename = data_filename
        self.map_filename = map_filename
        self.cols = cols

    def get_targets(self):
        """
        Parse data file and return targets and corresponding SMILES.
        """
        data = self.read_pcba_data(self.data_filename)
        id_map = read_pickle(self.map_filename)

        # get compound SMILES from map
        # indices are for data rows successfully mapped to SMILES
        smiles, indices = self.map_cids_to_smiles(data.PUBCHEM_CID, id_map)

        # get targets
        if self.cols is not None:
            cols = np.asarray(self.cols, dtype=int)
            targets = np.zeros((data.shape[0], len(cols)), dtype=float)
            for i, idx in enumerate(cols):
                targets[:, i] = data[data.columns[idx]]
        else:
            targets = np.asarray(data.PUBCHEM_ACTIVITY_OUTCOME == 'Active')
        targets = targets[indices]  # reduce targets to matched structures
        return smiles, targets

    def read_pcba_data(self, filename):
        """
        Read PCBA data.

        Parameters
        ----------
        filename : str
            PCBA data filename.
        """
        if filename.endswith('.gz'):
            with gzip.open(filename) as f:
                df = pd.read_csv(f)
        else:
            df = pd.read_csv(filename)
        df = df.drop_duplicates('PUBCHEM_CID')  # remove duplicate CIDs
        return df

    def map_cids_to_smiles(self, cids, id_map):
        """
        Look up SMILES for PubChem CIDs in a CID->SMILES map.

        Parameters
        ----------
        cids : array_like
            List of PubChem CIDs.
        id_map : dict
            CID->SMILES map.
        """
        smiles = []
        indices = []
        for i, cid in enumerate(cids):
            if np.isnan(cid):
                continue
            cid = int(cid)  # CIDs are often read in as floats
            name = 'CID{}'.format(cid)  # no bare IDs allowed in maps
            if name in id_map:
                smiles.append(id_map[name])
                indices.append(i)
        return np.asarray(smiles), np.asarray(indices)

    def get_column_names(self):
        """
        Get names of selected data columns.
        """
        if self.cols is None:
            return
        names = []
        for i in self.cols:
            names.append(self.read_pcba_data(self.data_filename).columns[i])
        return names
