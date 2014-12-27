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
    def __init__(self, features, targets):
        self.features = features
        self.targets = targets
        self.X = None
        self.y = None
        self.smiles = None
        self._get_dataset()

    def _get_dataset(self):
        """
        Get features, targets, and SMILES for this dataset.
        """
        # load features
        features = []
        feature_smiles = []
        for filename in self.features:
            data = read_pickle(filename)
            features.append(data['features'])
            feature_smiles.append(data['smiles'])
        features = np.ma.vstack(features)
        feature_smiles = np.concatenate(feature_smiles)

        # load targets
        data = read_pickle(self.targets)
        targets = data['targets']
        target_smiles = data['smiles']

        # match targets and features
        features_mask = np.in1d(feature_smiles, target_smiles)
        features_sort = np.argsort(feature_smiles[features_mask])
        targets_mask = np.in1d(target_smiles, feature_smiles)
        targets_sort = np.argsort(target_smiles[targets_mask])
        assert np.array_equal(feature_smiles[features_mask][features_sort],
                              target_smiles[targets_mask][targets_sort])
        self.X = features[features_mask][features_sort]
        self.y = targets[targets_mask][targets_sort]
        self.smiles = target_smiles[targets_mask][targets_sort]
