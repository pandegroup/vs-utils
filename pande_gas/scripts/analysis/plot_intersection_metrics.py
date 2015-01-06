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
    parser.add_argument('-i', '--inter', required=1, nargs='+',
                        help='Intersections.')
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
    return parser.parse_args(input_args)


def main(inter_filenames, scores_filename, output_filename, sim=False,
         actives=False, datasets=None):
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
    if datasets is not None:
        for dataset in datasets:
            m = re.search('^(.*?)-', os.path.basename(dataset))
            name, = m.groups()
            #name = name.replace('_', '-')  # fix for names
            data = h5_utils.load(dataset)
            targets[name] = data['y'][:]

    # get average intersection percentage for each dataset
    inter = {}
    for inter_filename in inter_filenames:
        m = re.search('^(.*?)-(.*?)-', os.path.basename(inter_filename))
        a, b = m.groups()
        if a == b:
            continue  # don't count self-intersection
        data = read_pickle(inter_filename)
        if 'sim' in data:
            assert data['sim'].shape == data['inter'].shape
        if a not in inter:
            inter[a] = []
        if actives:
            assert data['inter'].shape == targets[a].shape
            sel = np.where(targets[a])[0]
            if sim:
                inter[a].append(np.mean(data['sim'][sel]))
            else:
                inter[a].append(
                    np.count_nonzero(data['inter'][sel]) / np.count_nonzero(targets[a]))
        else:
            if sim:
                inter[a].append(np.mean(data['sim']))
            else:
                inter[a].append(np.count_nonzero(data['inter']) / data['inter'].size)

    for key in inter:
        print key, len(inter[key])
        #inter[key] = np.mean(inter[key])

    # get scores
    df = pd.read_table(scores_filename)
    scores = {}
    for name, score in zip(df.columns[1:], df.values[10][1:]):
        if name.startswith('PCBA'):
            name = name.split('PCBA-AID')[-1]
        elif name.startswith('MUV'):
            name = name.split('MUV-')[-1]
        elif name.startswith('TOX'):
            name = name.split('-')
            name.pop()
            name.pop(0)
            name = '_'.join(name)
        elif name.startswith('DUDE'):
            name = name.split('DUDE-')[-1]
        scores[name] = score

    assert np.all(np.in1d(inter.keys(), scores.keys()))

    # plot
    x, y, x_err = [], [], []
    for key in inter.keys():
        x.append(np.mean(inter[key]))
        x_err.append(np.std(inter[key]))
        y.append(scores[key])
    fig = pp.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y)
    ax.errorbar(x, y, xerr=x_err, linestyle='None')
    if sim:
        ax.set_xlabel('Mean Max Tanimoto Similarity')
    else:
        ax.set_xlabel('Mean Intersection')
    ax.set_ylabel('Mean AUC')
    fig.savefig(output_filename, dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    args = get_args()
    main(args.inter, args.scores, args.output, args.sim, args.actives, args.datasets)
