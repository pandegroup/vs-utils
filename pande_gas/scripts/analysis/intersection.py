"""
Look at the intersection between datasets.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import cPickle
import gzip
import matplotlib.pyplot as pp
import numpy as np
from scipy.misc import comb
import sys


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
    return parser.parse_args(input_args)


def main(input_filenames, output_prefix):
    """
    Look at the intersection between datasets.

    Parameters
    ----------
    input_filenames : list
        Input target filenames.
    output_prefix : str
        Prefix for output files.
    """
    n = len(input_filenames)
    smiles = []
    inter = np.ma.zeros((n, n), dtype=int)
    total = np.zeros(n, dtype=int)
    for filename in input_filenames:
        print filename
        with gzip.open(filename) as f:
            data = cPickle.load(f)
            smiles.append(data['smiles'])

    print 'Calculating {} intersections'.format(comb(n, 2, exact=True)),
    for i in xrange(n):
        total[i] = len(smiles[i])
        for j in xrange(n):
            if i == j:
                inter[i, j] = total[i]
                continue
            if i < j:
                continue
            inter[i, j] = np.intersect1d(smiles[i], smiles[j]).size
            inter[j, i] = inter[i, j]
        sys.stdout.write('.')
        sys.stdout.flush()
    print

    with gzip.open('{}-inter.pkl.gz'.format(output_prefix), 'wb') as f:
        cPickle.dump({'filenames': input_filenames, 'inter': inter,
                      'total': total}, f, cPickle.HIGHEST_PROTOCOL)

    # mask self-intersections
    inter.mask = np.zeros_like(inter)
    inter.mask[np.diag_indices_from(inter)] = True

    # calculate fractional intersections (breaks symmetry)
    inter = np.true_divide(inter, total)

    # plot
    fig = pp.figure()
    ax = fig.add_subplot(111)
    im = ax.imshow(inter, interpolation='None', origin='lower')
    fig.colorbar(im)
    pp.show()

if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.prefix)
