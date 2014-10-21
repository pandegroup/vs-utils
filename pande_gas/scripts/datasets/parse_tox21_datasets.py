#!/usr/bin/env python
"""
Get Tox21 challenge datasets.

Some SMILES strings are represented by multiple compounds, with some overlap of
assays. These compounds need to be condensed and their assay outcomes need to
be reconciled.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "3-clause BSD"

import argparse
import cPickle
import gzip
import numpy as np
import os
import warnings

from pande_gas.utils.target_utils import Tox21Parser


def get_args(input_args=None):
    """
    Get command-line arguments.

    Parameters
    ----------
    input_args : list, optional
        Input arguments. If not provided, defaults to sys.argv[1:].
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        help='Input filename.')
    parser.add_argument('--merge', choices=['max', 'min', 'majority_neg',
                                            'majority_pos'], default='max',
                        help='Target merge strategy.')
    parser.add_argument('-d', '--dir', default='.',
                        help='Directory in which to write target files.')
    return parser.parse_args(input_args)


def main(filename, merge_strategy, directory='.'):
    """
    Get Tox21 chellenge datasets.

    Parameters
    ----------
    filename : str
        Data filename.
    merge_strategy : str, optional (default 'max')
        Strategy to use when merging targets for duplicated molecules. Choose
        from 'max' (active if active in any assay), 'min' (inactive if inactive
        in any assay), 'majority_pos' (majority vote with ties assigned
        active), or 'majority_neg' (majority vote with ties assigned inactive).
    directory : str, optional (default '.')
        Directory in which to write target files.
    """
    parser = Tox21Parser(filename, merge_strategy=merge_strategy)
    data = parser.get_targets()

    # save individual datasets
    for dataset in data:
        if not len(data[dataset]):
            warnings.warn('Dataset "{}" is empty.'.format(dataset))
            continue
        targets = data[dataset]['targets']
        pos = np.count_nonzero(targets == 1)
        neg = np.count_nonzero(targets == 0)
        assert pos + neg == targets.size
        print '{}\t{}\t{}'.format(dataset, pos, neg)
        filename = os.path.join(directory, '{}-classes.pkl.gz'.format(dataset))
        with gzip.open(filename, 'wb') as f:
            cPickle.dump(data[dataset], f, cPickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
    args = get_args()
    main(args.input, args.merge, args.dir)
