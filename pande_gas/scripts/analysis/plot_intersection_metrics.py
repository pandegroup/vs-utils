"""
Plot intersection metrics.
"""
from __future__ import division

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pp
import numpy as np
import os
import re
from scipy.stats import linregress
import seaborn as sns
sns.set(style='whitegrid')
sns.set_palette('colorblind')
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
    parser.add_argument('--no-dude', action='store_true')
    parser.add_argument('--log-odds', action='store_true')
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


def plot_cor(inter, sizes, scores, datasets, output_filename, metric, no_dude=False, log_odds=False):
    # sanity checks
    for key in inter:
        #print key, len(inter[key])
        assert len(inter[key]) == sizes[key]
    assert np.all(np.in1d(inter.keys(), scores.keys()))

    # get x and y
    x, y, x_err = [], [], []
    names = []
    all_cor = []
    for key in inter.keys():
        if no_dude and key in datasets['DUDE']:
            continue
        names.append(key)
        x.append(np.mean(inter[key]))
        x_err.append(np.std(inter[key]))
        y.append(scores[key])
        if key in datasets['MUV']:
          all_cor.append(inter[key])
    x = np.asarray(x)
    x_err = np.asarray(x_err)
    y = np.asarray(y)
    names = np.asarray(names)
    print "DATASETS:", len(x)
    all_cor = np.concatenate(all_cor)
    fig = pp.figure()
    ax = fig.add_subplot(111)
    ax.hist(all_cor, bins=np.arange(30), normed=1)
    fig.savefig(output_filename + 'hist.png', dpi=300)

    with open('results.txt', 'wb') as f:
      for i, name in enumerate(names):
        if name.startswith('aid'):
            name = 'muv-' + name
        elif name.startswith('SR') or name.startswith('NR'):
            name = 'tox-' + name.replace('_', '-')
        elif name.startswith(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')):
            name = 'pcba-aid' + name
        else:
            name = 'dude-' + name
        f.write('{}\t{}\t{}\t{}\n'.format(name, x[i], x_err[i], y[i]))

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
        if no_dude and key == 'DUDE':
            print 'NOT PLOTTING DUDE'
            continue
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
    #ax.plot([-5, 25], [-5*m + b, 25*m + b])
    from sklearn.metrics import r2_score
    #print 'R2:', r2_score(x, m*x+b)
    ax.set_xlim(0, None)
    if metric == 'aor':
      ax.set_xlim(0, 20)
    if metric == 'cor':
        ax.set_xlabel(r'Compound Occurrence Rate (COR$_{i, \alpha}$)')
    elif metric == 'aor':
        ax.set_xlabel(r'Active Occurrence Rate (AOR$_{i, \alpha}$)', fontsize=14)
    elif metric == 'sim':
        ax.set_xlabel('Mean-Max Tanimoto Similarity')
    else:
        raise NotImplementedError(metric)
    if log_odds:
        ax.set_ylabel(r'$\Delta$ Log-Odds-Mean-AUC', fontsize=14)
    else:
        ax.set_ylabel(r'$\Delta$ Mean-AUC')
    pp.legend(loc=0, fontsize=14)
    fig.savefig(output_filename, dpi=300, bbox_inches='tight',
                transparent=True)


def plot_heatmap(inter_pairwise, datasets, output_filename):
    groups = OrderedDict()
    groups['PCBA'] = 'PCBA'
    groups['MUV'] = 'MUV'
    groups['Tox21'] = 'TOX'
    groups['DUD-E'] = 'DUDE'
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
            m[i, j] = inter_pairwise[a][b]
    assert np.count_nonzero(m.mask) == 0
    cmap = sns.cubehelix_palette(8, start=0.5, rot=-0.75, as_cmap=True)
    cmap = None
    cmap = 'RdBu_r'
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
              color='black')
    ax.hlines([259-128, 259-145, 259-157], 0, 259, linestyles='solid',
              linewidth=0.1, color='black')

    ax.text(50, 261, 'PCBA', fontdict={'fontsize': 6})
    ax.text(131, 261, 'MUV', fontdict={'fontsize': 6})
    ax.text(145, 261, 'Tox21', fontdict={'fontsize': 6})
    ax.text(195, 261, 'DUD-E', fontdict={'fontsize': 6})

    ax.text(-5, 205, 'PCBA', fontdict={'fontsize': 6, 'rotation': 90})
    ax.text(-5, 124, 'MUV', fontdict={'fontsize': 6, 'rotation': 90})
    ax.text(-5, 110, 'Tox21', fontdict={'fontsize': 6, 'rotation': 90})
    ax.text(-5, 60, 'DUD-E', fontdict={'fontsize': 6, 'rotation': 90})

    fig.savefig(output_filename, dpi=300, bbox_inches='tight',
                transparent=True)
    #import IPython
    #IPython.embed()


def main(inter_filenames, scores_filename, output_filename,
         target_filenames=None, no_dude=False, log_odds=False):
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
    scores, datasets = get_scores(scores_filename, log_odds=log_odds)

    if len(inter_filenames) == 1:
        (inter, inter_pairwise, active_inter,
            active_inter_pairwise) = read_pickle(inter_filenames[0])
    else:
        inter = {}
        inter_pairwise = {}
        active_inter = {}
        active_inter_pairwise = {}
        active_sim = {}
        active_sim_pairwise = {}
        for inter_filename in inter_filenames:
            m = re.search('^(.*?)-(.*?)-', os.path.basename(inter_filename))
            a, b = m.groups()
            data = read_pickle(inter_filename)
            a_inter = np.asarray(data['inter'], dtype=int)
            a_active_inter = np.asarray(data['active_inter'], dtype=int)
            a_active_sim = np.asarray(data['active_tanimoto'], dtype=float)

            # sanity checks
            assert active_sizes[a] == np.count_nonzero(targets[a])
            assert active_sizes[a] == a_active_inter.size
            assert sizes[a] == a_inter.size

            # get metric
            if a != b:  # don't count self-intersections
                if a not in inter:
                    inter[a] = np.zeros_like(a_inter)
                inter[a] += a_inter

                # "active" COR needs to measure the actives that are *active*
                # in other datasets
                # this requires new intersections, I think
                if a not in active_inter:
                    active_inter[a] = np.zeros_like(a_active_inter)
                active_inter[a] += a_active_inter

                if a not in active_sim:
                    active_sim[a] = []
                active_sim[a].append(a_active_sim)

            # get pairwise metric
            # fraction of dataset A in dataset B
            if a not in inter_pairwise:
                inter_pairwise[a] = {}
            assert b not in inter_pairwise[a]
            inter_pairwise[a][b] = np.true_divide(
                np.count_nonzero(a_inter), a_inter.size)

            if a not in active_inter_pairwise:
                active_inter_pairwise[a] = {}
            assert b not in active_inter_pairwise[a]
            active_inter_pairwise[a][b] = np.true_divide(
                np.count_nonzero(a_active_inter), a_active_inter.size)

            if a not in active_sim_pairwise:
                active_sim_pairwise[a] = {}
            assert b not in active_sim_pairwise[a]
            active_sim_pairwise[a][b] = np.mean(a_active_sim)

        # correct active_sim to be n_compounds x n_assays
        for key in active_sim.iterkeys():
            active_sim[key] = np.asarray(active_sim[key], dtype=float).T

        write_pickle(
            [inter, inter_pairwise, active_inter, active_inter_pairwise],
            'data.pkl.gz')

    plot_cor(inter, sizes, scores, datasets, 'cor.png', 'cor', no_dude, log_odds)
    plot_cor(active_inter, active_sizes, scores, datasets, 'cor-actives.png', 'aor', no_dude, log_odds)
    #plot_cor(active_sim, active_sizes, scores, datasets, 'cor-sim.png', 'sim', no_dude, log_odds)

    plot_heatmap(inter_pairwise, datasets, 'heatmap.png')
    plot_heatmap(active_inter_pairwise, datasets, 'heatmap-actives.png')
    #plot_heatmap(active_sim_pairwise, datasets, 'heatmap-sim.png')

if __name__ == '__main__':
    args = get_args()
    if args.inter:
        inter = args.inter
    else:
        inter = []
        with open(args.file) as f:
            for line in f:
                inter.append(line.strip())
    main(inter, args.scores, args.output, args.targets, args.no_dude, args.log_odds)
