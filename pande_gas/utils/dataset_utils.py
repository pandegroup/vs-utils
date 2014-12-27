"""
Dataset preparation utilities.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import gzip
import numpy as np

from rdkit import Chem

from . import read_pickle, SmilesGenerator


class MoleculeDatabase(object):
    """
    Molecule database.

    Molecules are keyed by SMILES.

    Parameters
    ----------
    kwargs : dict, optional
        Keyword arguments for SmilesMap.
    """
    def __init__(self, **kwargs):
        self.engine = SmilesGenerator(**kwargs)
        self.smiles = set()

    def __len__(self):
        return len(self.smiles)

    def __iter__(self):
        for smiles in self.smiles:
            yield smiles

    def __contains__(self, item):
        return item in self.smiles

    def add_mol(self, mol):
        """
        Add a molecule to the database.

        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        self.smiles.add(self.engine.get_smiles(mol))

    def load(self, filename):
        """
        Load an existing database.

        Parameters
        ----------
        filename : str
            Existing database filename.
        """
        if filename.endswith('.gz'):
            f = gzip.open(filename)
        else:
            f = open(filename)
        for line in f:
            smiles = line.strip()
            mol = Chem.MolFromSmiles(smiles)  # sanity check
            if mol is None:
                raise ValueError(
                    'Database is unreadable: "{}".'.format(smiles))
            self.smiles.add(smiles)
        f.close()

    def save(self, filename):
        """
        Save the database to disk.

        Parameters
        ----------
        filename : str
            Filename.
        """
        if filename.endswith('.gz'):
            f = gzip.open(filename, 'wb')
        else:
            f = open(filename, 'wb')
        for smiles in self.smiles:
            f.write('{}\n'.format(smiles))
        f.close()


class Dataset(object):
    """
    Extract features for a specific dataset.

    Parameters
    ----------
    features : list
        Feature filenames.
    targets : str
        Target filename.
    """
    def __init__(self, features, feature_smiles, targets, target_smiles):
        self.X, self.y, self.smiles = self.get_dataset(
            features, feature_smiles, targets, target_smiles)

    def get_dataset(self, features, feature_smiles, targets, target_smiles):
        """
        Get features, targets, and SMILES for a dataset.

        Parameters
        ----------
        features : array_like
            Features.
        feature_smiles : array_like
            SMILES for features.
        targets : array_like
            Targets.
        target_smiles : array_like
            SMILES for targets.

        Returns
        -------
        X : array_like
            Features matching targets.
        y : array_like
            Targets.
        smiles : array_like
            SMILES matching targets.
        """
        # match targets and features
        # there may be duplicate target SMILES, so use searchsorted
        features_mask = np.in1d(feature_smiles, target_smiles)
        features_sort = np.argsort(feature_smiles[features_mask])
        targets_mask = np.in1d(target_smiles, feature_smiles)
        targets_sort = np.argsort(target_smiles[targets_mask])
        features_sel = np.searchsorted(
            feature_smiles[features_mask][features_sort],
            target_smiles[targets_mask][targets_sort])
        assert np.array_equal(
            feature_smiles[features_mask][features_sort][features_sel],
            target_smiles[targets_mask][targets_sort])
        X = features[features_mask][features_sort][features_sel]
        y = targets[targets_mask][targets_sort]
        smiles = target_smiles[targets_mask][targets_sort]
        return X, y, smiles
