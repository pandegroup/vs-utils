"""
Calculate similarity between two datasets using feature Tanimotos.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import numpy as np

from pande_gas.utils import h5_utils


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
    return parser.parse_args(input_args)


def main(ref_filename, fit_filename, output_filename):
    """
    Calculate similarity between two datasets.

    Parameters
    ----------
    ref_filename : str
        Reference dataset.
    fit_filename : str
        Fit dataset.
    output_filename : str
        Output filename.
    """
    ref = h5_utils.load(ref_filename)
    fit = h5_utils.load(fit_filename)
    sim = tanimoto(ref['X'][:], fit['X'][:])
    h5_utils.dump(sim, output_filename)


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
    ref_self_overlap = np.diag(np.dot(ref, ref.T))
    fit_self_overlap = np.diag(np.dot(fit, fit.T))
    ref_fit_overlap = np.dot(ref, fit.T)
    return np.true_divide(
        ref_fit_overlap, ref_self_overlap + fit_self_overlap - ref_fit_overlap)

if __name__ == '__main__':
    args = get_args()
    main(args.ref, args.fit, args.output)
