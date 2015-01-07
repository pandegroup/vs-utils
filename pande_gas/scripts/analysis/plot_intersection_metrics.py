"""
Plot intersection metrics.
"""
from __future__ import division

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pp
import numpy as np
import os
import pandas as pd
import re
from scipy.stats import linregress
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

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
    parser.add_argument('--actives', action='store_true',
                        help='Only use actives in metrics.')
    parser.add_argument('-d', '--datasets', nargs='+',
                        help='Datasets containing labels.')
    parser.add_argument('-o', '--output', required=1,
                        help='Output filename.')
    parser.add_argument('--sim', action='store_true',
                        help='Calculate similarity metrics.')
    parser.add_argument('--union', action='store_true',
                        help='Compare vs. union of all other datasets.')
    return parser.parse_args(input_args)


def get_weighted_mean_and_std(x, w):
    """
    Weighted mean and standard deviation.

    See http://stats.stackexchange.com/questions/6534.

    Parameters
    ----------
    x : array_like
        Observations.
    w : array_like
        Weights.
    """
    mean = np.true_divide(np.sum(np.multiply(w, x)), np.sum(w))
    nonzero = np.count_nonzero(w)
    stdev = np.sqrt(np.true_divide(
        np.sum(np.multiply(w, np.square(x - mean))),
        np.true_divide(nonzero - 1, nonzero) * np.sum(w)))
    return mean, stdev


def main(inter_filenames, scores_filename, output_filename, sim=False,
         actives=False, datasets=None, union=False):
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
    targets = {}
    sizes = {}
    active_sizes = {}
    if datasets is not None:
        for dataset in datasets:
            m = re.search('^(.*?)-', os.path.basename(dataset))
            name, = m.groups()
            #name = name.replace('_', '-')  # fix for names
            data = h5_utils.load(dataset)
            targets[name] = data['y'][:]
            sizes[name] = data['X'].shape[0]
            active_sizes[name] = np.count_nonzero(data['y'])

    # get average intersection percentage for each dataset
    inter = {}
    for inter_filename in inter_filenames:
        m = re.search('^(.*?)-(.*?)-', os.path.basename(inter_filename))
        if m is None:
            import IPython
            IPython.embed()
            sys.exit()
        a, b = m.groups()
        if a == b:
            continue  # don't count self-intersection
        data = read_pickle(inter_filename)
        if 'sim' in data:
            assert data['sim'].shape == data['inter'].shape
        if a not in inter:
            if union:
                inter[a] = np.zeros_like(data['inter'], dtype=float)
            else:
                inter[a] = {}
        assert b not in inter[a]

        # get metric
        if actives:
            assert data['inter'].shape == targets[a].shape
            a_sel = np.where(targets[a])[0]
            b_sel = np.where(targets[b])[0]
            if sim:
                if union:
                    inter[a] = np.maximum(inter[a], data['sim'])
                else:
                    inter[a][b] = np.mean(data['sim'][sel])
            else:
                if union:
                    inter[a] = np.add(inter[a], data['inter'])
                else:
                    inter[a][b] = np.true_divide(
                        np.count_nonzero(data['inter'][sel]),
                        active_sizes[a])
        else:
            if sim:
                if union:
                    inter[a] = np.maximum(inter[a], data['sim'])
                else:
                    inter[a][b] = np.mean(data['sim'])
            else:
                if union:
                    inter[a] = np.add(inter[a], data['inter'])
                else:
                    inter[a][b] = np.true_divide(
                        np.count_nonzero(data['inter']),
                        data['inter'].size)
        if a in sizes:
            assert active_sizes[a] == np.count_nonzero(targets[a])
            assert sizes[a] == data['inter'].size

    for key in inter:
        print key, len(inter[key])
        if union:
            #if actives:
            #    assert len(inter[key]) == active_sizes[key]
            #else:
            #    assert len(inter[key]) == sizes[key]
            assert len(inter[key]) == sizes[key]
        else:
            assert len(inter[key]) == 258, (key, len(inter[key]))
        #inter[key] = np.mean(inter[key])

    # get scores
    df = pd.read_table(scores_filename)
    scores = {}
    ref_idx = 1
    new_idx = 27
    datasets = {'PCBA': [], 'MUV': [], 'TOX': [], 'DUDE': []}
    print df.values[ref_idx][0], df.values[new_idx][0]  # print scores
    for name, ref_score, new_score in zip(
            df.columns[1:], df.values[ref_idx][1:], df.values[new_idx][1:]):
        score = new_score - ref_score
        if name.startswith('PCBA'):
            name = name.split('PCBA-AID')[-1]
            datasets['PCBA'].append(name)
        elif name.startswith('MUV'):
            name = name.split('MUV-')[-1]
            datasets['MUV'].append(name)
        elif name.startswith('TOX'):
            name = name.split('-')
            name.pop()
            name.pop(0)
            name = '_'.join(name)
            datasets['TOX'].append(name)
        elif name.startswith('DUDE'):
            name = name.split('DUDE-')[-1]
            datasets['DUDE'].append(name)
        else:
            raise ValueError(name)
        scores[name] = score
    total = 0
    for key in datasets:
        total += len(datasets[key])
    assert total == 259, total

    assert np.all(np.in1d(inter.keys(), scores.keys()))

    # plot
    x, y, x_err = [], [], []
    names = []
    for key in inter.keys():
        names.append(key)
        if union:

            # fix active selection
            if actives:
                sel = np.where(targets[key])[0]
                inter[key] = inter[key][sel]

            # get metrics
            if sim:
                x.append(np.mean(inter[key]))
            else:
                x.append(np.count_nonzero(inter[key]) / inter[key].size)
        else:
            values, weights = [], []
            for other in inter[key].keys():
                values.append(inter[key][other])
                if actives:
                    weights.append(active_sizes[other])
                else:
                    weights.append(sizes[other])
            mean, stdev = get_weighted_mean_and_std(values, weights)
            x.append(mean)
            x_err.append(stdev)
        y.append(scores[key])
    x = np.asarray(x)
    x_err = np.asarray(x_err)
    y = np.asarray(y)
    names = np.asarray(names)

    # linear fit
    #lr = LinearRegression()
    #lr.fit(x, y)
    #r2 = r2_score(y, lr.predict(x))
    #print "R2", r2, r2_score(x, y)
    m, b, r, p, err = linregress(x, y)
    print m, b, r, p, err
    print "R2", r ** 2

    fig = pp.figure()
    ax = fig.add_subplot(111)
    #ax.plot([0, 1], [b, m + b], color='k')
    for key in datasets:
        if key == 'DUDE':
            print 'NOT PLOTTING DUDE'
            continue
        sel = []
        for name in datasets[key]:
            if name not in names:
                import IPython
                IPython.embed()
                sys.exit()
            try:
                idx = np.where(names == name)[0][0]
            except IndexError:
                import IPython
                IPython.embed()
                sys.exit()
            sel.append(idx)
        sel = np.asarray(sel, dtype=int)
        assert sel.size == len(datasets[key])
        ax.plot(x[sel], y[sel], 'o', label=key)
        if len(x_err):
            ax.errorbar(x[sel], y[sel], xerr=x_err[sel], linestyle='None')
    #ax.scatter(x, y)
    #if len(x_err):
    #    ax.errorbar(x, y, xerr=x_err, linestyle='None')
    if sim:
        ax.set_xlabel('Mean Max Tanimoto Similarity')
    else:
        ax.set_xlabel('Mean Intersection')
    ax.set_ylabel(r'$\Delta$ Mean AUC')
    pp.legend()
    fig.savefig(output_filename, dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    args = get_args()
    if args.inter:
        inter = args.inter
    else:
        inter = []
        with open(args.file) as f:
            for line in f:
                inter.append(line.strip())
    main(inter, args.scores, args.output, args.sim, args.actives,
         args.datasets, args.union)
