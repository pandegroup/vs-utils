"""
Calculate similarity between two datasets using feature Tanimotos.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
from joblib import delayed, Parallel
import numpy as np

from pande_gas.utils import h5_utils, write_pickle


def get_args(input_args=None):
    """
    Get command-line arguments.

    Parameters
    ----------
    input_args : list, optional
        Input arguments. If not provided, defaults to sys.argv[1:].
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', required=1,
                        help='Reference dataset.')
    parser.add_argument('-f', '--fit', required=1,
                        help='Fit dataset.')
    parser.add_argument('-o', '--output', required=1,
                        help='Output filename.')
    parser.add_argument('-np', '--n-jobs', type=int, default=1,
                        help='Number of parallel jobs.')
    return parser.parse_args(input_args)


def main(ref_filename, fit_filename, output_filename, chunk_size=1000,
         n_jobs=1):
    """
    Calculate similarity between two datasets.

    Report the average maximum similarity of a compound in dataset A to all
    compounds in dataset B.

    Parameters
    ----------
    ref_filename : str
        Reference dataset.
    fit_filename : str
        Fit dataset.
    output_filename : str
        Output filename.
    n_jobs : int, optional (default 1)
    """
    ref_data = h5_utils.load(ref_filename)
    fit_data = h5_utils.load(fit_filename)
    ref = ref_data['X']#[:]
    fit = fit_data['X'][:]
    '''
    rval = Parallel(n_jobs=n_jobs, verbose=5)(
        delayed(get_max_tanimoto)(this_ref, fit)
        for this_ref in np.array_split(ref, ref.shape[0] / chunk_size))
    '''
    rval = []
    i = 0
    while i < ref.shape[0]:
        print i
        this_ref = ref[i:i+chunk_size]
        rval.append(get_max_tanimoto(this_ref, fit))
        print len(rval)
    rval = np.asarray(rval, dtype=float)
    write_pickle(rval, output_filename)


def get_max_tanimoto(ref, fit):
    """
    Get the maximum Tanimoto between a reference molecule and another dataset.

    Parameters
    ----------
    ref : array_like
        Features for reference molecule.
    fit : array_like
        Features for fit molecules.
    """
    assert ref.ndim == 2
    return np.amax(tanimoto(ref, fit), axis=1)


def tanimoto(ref, fit):
    """
    Calculate Tanimotos.

    Parameters
    ----------
    ref : array_like
        Reference dataset features.
    fit : array_like
        Fit dataset features.
    """
    a_overlap = np.sum(np.multiply(ref, ref), axis=1)
    b_overlap = np.sum(np.multiply(fit, fit), axis=1)
    b_overlap, a_overlap = np.meshgrid(b_overlap, a_overlap)
    ab_overlap = np.dot(ref, fit.T)
    return np.true_divide(
        ab_overlap, a_overlap + b_overlap - ab_overlap)

if __name__ == '__main__':
    args = get_args()
    main(args.ref, args.fit, args.output)
