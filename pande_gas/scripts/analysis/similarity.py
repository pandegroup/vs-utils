"""
Calculate similarity between two datasets using feature Tanimotos.

First check for overlap between datasets and only calculate Tanimotos for
compounds that do not overlap. Then we get both metrics at once and reduce
total computation.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
from joblib import delayed, Parallel
import numpy as np

from pande_gas.utils import h5_utils, vector_tanimoto, write_pickle


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
    parser.add_argument('--chunk-size', type=int, default=100,
                        help='Parallel chunk size.')
    parser.add_argument('--identity', action='store_true',
                        help='Only compute identity metrics.')
    parser.add_argument('-np', '--n-jobs', type=int, default=1,
                        help='Number of parallel jobs.')
    return parser.parse_args(input_args)


def main(ref_filename, fit_filename, output_filename, chunk_size=100,
         identity=False, n_jobs=1):
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
    chunk_size : int, optional (default 100)
        Chunk size for Tanimoto calculations.
    identity : bool, optional (default False)
        Only compute identity metrics.
    n_jobs : int, optional (default 1)
    """
    ref_data = h5_utils.load(ref_filename)
    fit_data = h5_utils.load(fit_filename)

    ref_inter = np.in1d(ref_data['smiles'], fit_data['smiles'])
    print 'Ref intersection: {}/{}'.format(
        np.count_nonzero(ref_inter), ref_inter.size)

    a = np.where(ref_data['y'][:])[0]
    b = np.where(fit_data['y'][:])[0]
    ref_active_inter = np.in1d(ref_data['smiles'][:][a],
                               fit_data['smiles'][:][b])

    # compute active--active similarity
    a_fp = ref_data['X'][:][a]
    b_fp = fit_data['X'][:][b]
    active_tan = vector_tanimoto(a_fp, b_fp)
    #a_so = np.sum(np.multiply(a_fp, a_fp), axis=1)
    #b_so = np.sum(np.multiply(b_fp, b_fp), axis=1)
    #tan = vector_tanimoto(a_fp, b_fp, a_so, b_so)

    if identity:
        write_pickle({'inter': ref_inter, 'active_inter': ref_active_inter, 
                      'active_tanimoto': active_tan}, output_filename)
        return

    sel = ~ref_inter
    ref = ref_data['X'][:]
    fit = fit_data['X'][:]
    if np.array_equal(np.unique(fit), [0, 1]):
        fit_overlap = np.sum(fit, axis=1)
    else:
        fit_overlap = np.sum(np.multiply(fit, fit), axis=1)

    n_todo = np.count_nonzero(sel)
    print 'Need to calculate similarity for {} compounds'.format(n_todo)
    if n_todo:
        rval = Parallel(n_jobs=n_jobs, verbose=5)(
            delayed(get_max_tanimoto)(this_ref, fit, fit_overlap=fit_overlap)
            for this_ref in np.array_split(
                ref[sel], n_todo / chunk_size))

        sim = np.zeros_like(ref_inter, dtype=float)
        sim[ref_inter] = 1.0
        sim[sel] = np.concatenate(rval)
    else:
        sim = np.ones_like(ref_inter, dtype=float)

    print 'Ref similarity: {}/{}'.format(np.mean(sim), np.median(sim))

    write_pickle({'sim': sim, 'inter': ref_inter}, output_filename)


def get_max_tanimoto(ref, fit, ref_overlap=None, fit_overlap=None):
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
    return np.amax(vector_tanimoto(ref, fit, ref_overlap, fit_overlap), axis=1)

if __name__ == '__main__':
    args = get_args()
    main(args.ref, args.fit, args.output, args.chunk_size, args.n_jobs)
