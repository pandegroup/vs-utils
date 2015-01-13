"""
Plot intersection metrics.
"""
from __future__ import division

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
from matplotlib import pyplot as pp
import numpy as np
import os
import re
from scipy.stats import linregress
import seaborn as sns  # import for style

from . import get_scores

from pande_gas.utils import h5_utils, read_pickle


def get_args(input_args=None):
    """
    Get command-line arguments.

    Parameters
    ----------
    input_args : list, optional
        Input arguments. If not provided, defaults to sys.argv[1:].
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inter', required=0, nargs='+',
                        help='Intersections.')
    parser.add_argument('-f', '--file', help='Input filename.')
    parser.add_argument('-s', '--scores', required=1,
                        help='Scores.')
    parser.add_argument('-t', '--targets', nargs='+', required=1,
                        help='Datasets containing targets.')
    parser.add_argument('-o', '--output', required=1,
                        help='Output filename.')
    return parser.parse_args(input_args)


def get_targets(filenames):
    """
    Get targets and dataset sizes.
    """
    targets = {}
    sizes = {}
    active_sizes = {}
    for filename in filenames:
        m = re.search('^(.*?)-', os.path.basename(filename))
        name, = m.groups()
        data = h5_utils.load(filename)
        targets[name] = data['y'][:]
        sizes[name] = data['X'].shape[0]
        active_sizes[name] = np.count_nonzero(data['y'])
    return targets, sizes, active_sizes


def main(inter_filenames, scores_filename, output_filename,
         target_filenames=None):
    """
    Plot intersection metrics.

    Parameters
    ----------
    inter_filenames : list
        Intersections.
    scores_filename : str
        Scores.
    output_filename : str
        Output filename.
    sim : bool, optional (default False)
        Calculate similarity metrics.
    actives : bool, optional (default False)
        Only use actives in metrics (reqires datasets).
    datasets : list, optional
        Datasets containing labels.
    """
    targets, sizes, active_sizes = get_targets(target_filenames)
    scores, datasets = get_scores(scores_filename)

    inter = {}
    for inter_filename in inter_filenames:
        m = re.search('^(.*?)-(.*?)-', os.path.basename(inter_filename))
        a, b = m.groups()
        if a == b:
            continue  # don't count self-intersections
        data = read_pickle(inter_filename)

        # sanity checks
        assert active_sizes[a] == np.count_nonzero(targets[a])
        assert sizes[a] == data['inter'].size

        # get metric
        if a not in inter:
            inter[a] = np.zeros_like(data['inter'], dtype=int)
        inter[a] += np.asarray(data['inter'], dtype=int)

    # sanity checks
    for key in inter:
        #print key, len(inter[key])
        assert len(inter[key]) == sizes[key]
    assert np.all(np.in1d(inter.keys(), scores.keys()))

    # get x and y
    x, y, x_err = [], [], []
    names = []
    for key in inter.keys():
        if key in datasets['DUDE']:
            continue
        names.append(key)
        x.append(np.mean(inter[key]))
        x_err.append(np.std(inter[key]))
        y.append(scores[key])
    x = np.asarray(x)
    x_err = np.asarray(x_err)
    y = np.asarray(y)
    names = np.asarray(names)
    print "DATASETS:", len(x)

    # statistics
    m, b, r, p, err = linregress(x, y)
    print m, b, r, p, err
    print r, r**2

    # plot
    fig = pp.figure()
    ax = fig.add_subplot(111)
    for key in ['PCBA', 'DUDE', 'MUV', 'TOX']:
        if key == 'DUDE':
            print 'NOT PLOTTING DUDE'
            continue
        sel = []
        for name in datasets[key]:
            idx = np.where(names == name)[0][0]
            sel.append(idx)
        sel = np.asarray(sel, dtype=int)
        assert sel.size == len(datasets[key])  # sanity check
        ax.plot(x[sel], y[sel], 'o', label=key)
        if len(x_err):
            color = 'gray'
            #if key == 'PCBA':
            #    color = 'blue'
            #elif key == 'MUV':
            #    color = 'green'
            #elif key == 'TOX':
            #    color = 'red'
            ax.errorbar(x[sel], y[sel], xerr=x_err[sel], linestyle='None',
                        color=color, elinewidth=0.5)
    #ax.plot([-20, 140], [-20*m + b, 140*m + b])
    ax.set_xlabel('Compound Occurrence Rate')
    ax.set_ylabel(r'$\Delta$ Mean AUC')
    pp.legend(loc=0)
    fig.savefig(output_filename, dpi=300, bbox_inches='tight',
                transparent=True)

if __name__ == '__main__':
    args = get_args()
    if args.inter:
        inter = args.inter
    else:
        inter = []
        with open(args.file) as f:
            for line in f:
                inter.append(line.strip())
    main(inter, args.scores, args.output, args.targets)
