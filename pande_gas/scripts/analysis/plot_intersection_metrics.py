"""
Plot intersection metrics.
"""
from __future__ import division

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
from collections import OrderedDict
from matplotlib import pyplot as pp
import numpy as np
import os
import re
from scipy.stats import linregress
import seaborn as sns
sns.set(style='whitegrid')
#sns.set(context='paper')

from pande_gas.scripts.analysis import get_scores
from pande_gas.utils import h5_utils, read_pickle, write_pickle


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

    if len(inter_filenames) == 1:
        inter, inter_pairwise = read_pickle(inter_filenames[0])
    else:
        inter = {}
        inter_pairwise = {}
        for inter_filename in inter_filenames:
            m = re.search('^(.*?)-(.*?)-', os.path.basename(inter_filename))
            a, b = m.groups()
            data = read_pickle(inter_filename)

            # sanity checks
            assert active_sizes[a] == np.count_nonzero(targets[a])
            assert sizes[a] == data['inter'].size

            # get metric
            if a != b:  # don't count self-intersections
                if a not in inter:
                    inter[a] = np.zeros_like(data['inter'], dtype=int)
                inter[a] += np.asarray(data['inter'], dtype=int)

            # get pairwise metric
            if a not in inter_pairwise:
                inter_pairwise[a] = {}
            assert b not in inter_pairwise[a]
            inter_pairwise[a][b] = np.count_nonzero(
                np.asarray(data['inter'], dtype=int))
        write_pickle([inter, inter_pairwise], 'data.pkl.gz')

    # sanity checks
    for key in inter:
        #print key, len(inter[key])
        assert len(inter[key]) == sizes[key]
    assert np.all(np.in1d(inter.keys(), scores.keys()))

    # get x and y
    x, y, x_err = [], [], []
    names = []
    for key in inter.keys():
        #if key in datasets['DUDE']:
        #    continue
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
    groups = OrderedDict()
    groups['PCBA'] = 'PCBA'
    groups['MUV'] = 'MUV'
    groups['Tox21'] = 'TOX'
    groups['DUD-E'] = 'DUDE'
    for group, key in groups.iteritems():
        #if key == 'DUDE':
        #    print 'NOT PLOTTING DUDE'
        #    continue
        sel = []
        for name in datasets[key]:
            try:
                idx = np.where(names == name)[0][0]
            except IndexError:
                import IPython
                IPython.embed()
                sys.exit()
            sel.append(idx)
        sel = np.asarray(sel, dtype=int)
        assert sel.size == len(datasets[key])  # sanity check
        ax.plot(x[sel], y[sel], 'o', label=group)
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
    ax.set_xlim(0, None)
    ax.set_xlabel(r'Compound Occurrence Rate (COR$_{i, \alpha}$)')
    ax.set_ylabel(r'$\Delta$ Mean AUC')
    pp.legend(loc=0)
    fig.savefig(output_filename, dpi=300, bbox_inches='tight',
                transparent=True)

    # plot full pairwise heatmap
    # the color bar will show the fraction of dataset A in dataset B
    fig = pp.figure()
    ax = fig.add_subplot(111)
    m = np.ma.masked_all(
        (len(inter_pairwise), len(inter_pairwise)), dtype=float)
    names = []
    for group, key in groups.iteritems():
        for name in datasets[key]:
            names.append(name)
    names = np.asarray(names)
    print names
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            m[i, j] = np.true_divide(inter_pairwise[a][b], sizes[a])
    assert np.count_nonzero(m.mask) == 0
    #cmap = sns.cubehelix_palette(8, start=0.5, rot=-0.75, as_cmap=True)
    cmap = None
    '''
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import linkage
    s = np.zeros_like(m, dtype=float)
    for group, key in groups.iteritems():
        for name in datasets[key]:
            i = np.where(names == name)[0][0]
            for j, other in enumerate(names):
                if other in datasets[key]:
                    s[i, j] = 0
                else:
                    s[i, j] = 1
    link = linkage(squareform(s))
    sns.clustermap(m, square=True, xticklabels=False, yticklabels=False,
                   linewidths=0, cmap=cmap, row_linkage=link, col_linkage=link)
    '''
    sns.heatmap(m, square=True, xticklabels=False, yticklabels=False,
                linewidths=0, cmap=cmap)
    ax.vlines([128, 145, 157], 0, 259, linestyles='solid', linewidth=0.1,
              color='red')
    ax.hlines([259-128, 259-145, 259-157], 0, 259, linestyles='solid',
              linewidth=0.1, color='red')

    ax.text(50, 261, 'PCBA', fontdict={'fontsize': 6})
    ax.text(131, 261, 'MUV', fontdict={'fontsize': 6})
    ax.text(145, 261, 'Tox21', fontdict={'fontsize': 6})
    ax.text(195, 261, 'DUD-E', fontdict={'fontsize': 6})

    ax.text(-5, 205, 'PCBA', fontdict={'fontsize': 6, 'rotation': 90})
    ax.text(-5, 124, 'MUV', fontdict={'fontsize': 6, 'rotation': 90})
    ax.text(-5, 110, 'Tox21', fontdict={'fontsize': 6, 'rotation': 90})
    ax.text(-5, 60, 'DUD-E', fontdict={'fontsize': 6, 'rotation': 90})

    fig.savefig('heatmap.png', dpi=300, bbox_inches='tight', transparent=True)
    #import IPython
    #IPython.embed()

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
