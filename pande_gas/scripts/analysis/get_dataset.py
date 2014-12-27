"""
Get features and targets for one or more datasets.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import numpy as np
import os

from pande_gas.utils import read_pickle, write_pickle
from pande_gas.utils.dataset_utils import Dataset


def get_args(input_args=None):
    """
    Get command-line arguments.

    Parameters
    ----------
    input_args : list, optional
        Input arguments. If not provided, defaults to sys.argv[1:].
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--features', required=1, nargs='+',
                        help='Features.')
    parser.add_argument('-t', '--targets', required=1, nargs='+',
                        help='Targets.')
    parser.add_argument('-s', '--suffix', required=1,
                        help='Suffix for output filenames.')
    return parser.parse_args(input_args)


def main(feature_filenames, target_filenames, suffix):
    """
    Get features and targets for one or more datasets.

    Parameters
    ----------
    feature_filenames : list
        Feature filename(s).
    target_filenames : list
        Targets filename(s).
    suffix : str
        Suffix for output filenames.
    """
    # load features
    features = []
    feature_smiles = []
    for feature_filename in feature_filenames:
        data = read_pickle(feature_filename)
        features.append(data['features'])
        feature_smiles.append(data['smiles'])
    features = np.ma.vstack(features)
    feature_smiles = np.concatenate(feature_smiles)

    for target_filename in target_filenames:
        prefix = os.path.basename(target_filename).split('-')[0]

        # load targets
        data = read_pickle(target_filename)
        targets = data['targets']
        target_smiles = data['smiles']

        dataset = Dataset(features, feature_smiles, targets, target_smiles)
        print '{}/{} target records extracted'.format(
            dataset.X.shape[0], targets.shape[0])
        write_pickle(dataset, '{}-{}.pkl.gz'.format(prefix, suffix))

if __name__ == '__main__':
    args = get_args()
    main(args.features, args.targets, args.suffix)
