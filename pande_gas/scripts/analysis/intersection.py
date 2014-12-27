"""
Look at the intersection between datasets.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import cPickle
import gzip
from joblib import delayed, Parallel
import matplotlib.pyplot as pp
import numpy as np
from scipy.misc import comb


def parse_args(input_args=None):
    """
    Parse command-line arguments.

    Parameters
    ----------
    input_args : list, optional
        Input arguments. If not provided, defaults to sys.argv[1:].
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=1, nargs='+',
                        help='Input target file(s).')
    parser.add_argument('-p', '--prefix', required=1,
                        help='Prefix for output files.')
    parser.add_argument('-np', '--n-jobs', type=int, default=1,
                        help='Number of jobs.')
    return parser.parse_args(input_args)


def main(input_filenames, output_prefix, n_jobs=1):
    """
    Look at the intersection between datasets.

    Parameters
    ----------
    input_filenames : list
        Input target filenames.
    output_prefix : str
        Prefix for output files.
    n_jobs : int, optional (default 1)
        Number of jobs.
    """
    n = len(input_filenames)
    smiles = []
    inter = np.ma.zeros((n, n), dtype=int)
    total = np.zeros(n, dtype=int)
    for i, filename in enumerate(input_filenames):
        print filename
        if filename.endswith('.gz'):
            f = gzip.open(filename)
        else:
            f = open(filename)
        data = cPickle.load(f)
        f.close()
        smiles.append(data['smiles'])
        total[i] = len(data['smiles'])

    print 'Calculating {} intersections'.format(comb(n, 2, exact=True))

    # calculate upper triangular portion (minus diagonals)
    for i in xrange(n):
        inter[i, i+1:] = Parallel(n_jobs=n_jobs, verbose=5)(
            delayed(get_intersection_size)(smiles[i], this_smiles)
            for this_smiles in smiles[i+1:])

    # fill in remaining elements
    inter += inter.T
    inter[np.diag_indices_from(inter)] = total

    with gzip.open('{}-inter.pkl.gz'.format(output_prefix), 'wb') as f:
        cPickle.dump({'filenames': input_filenames, 'inter': inter,
                      'total': total}, f, cPickle.HIGHEST_PROTOCOL)

    # mask self-intersections
    # inter.mask = np.zeros_like(inter)
    # inter.mask[np.diag_indices_from(inter)] = True

    # calculate fractional intersections (breaks symmetry)
    inter = np.true_divide(inter, total)

    # plot
    fig = pp.figure()
    ax = fig.add_subplot(111)
    im = ax.imshow(inter, interpolation='None', origin='lower')
    fig.colorbar(im)
    pp.show()


def get_intersection_size(a, b):
    """
    Get the intersection size between two data sets.

    Parameters
    ----------
    a, b : array_like
        Datasets to compare.
    """
    return np.intersect1d(a, b).size

if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.prefix, args.n_jobs)
