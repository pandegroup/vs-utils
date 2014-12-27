"""
Get features and targets for a dataset.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse

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
    parser.add_argument('-t', '--targets', required=1,
                        help='Targets.')
    parser.add_argument('-o', '--output', required=1,
                        help='Output filename.')
    return parser.parse_args(input_args)


def main(features, targets, output):
    """
    Get features and targets for a dataset.

    Parameters
    ----------
    features : list
        Feature filename(s).
    targets : str
        Targets filename.
    output : str
        Output filename.
    """
    dataset = Dataset(features, targets)
    target_data = read_pickle(targets)
    print '{}/{} target records extracted'.format(
        dataset.X.shape[0], target_data['targets'].shape[0])
    write_pickle(dataset, output)

if __name__ == '__main__':
    main(**vars(get_args()))
