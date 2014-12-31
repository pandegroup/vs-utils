"""
Plot intersection metrics.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
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
    return parser.parse_args(input_args)


def main(inter_filenames, scores_filename, output_filename, actives=False,
         datasets=None):
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
    actives : bool, optional (default False)
        Only use actives in metrics (reqires datasets).
    datasets : list, optional
        Datasets containing labels.
    """
    targets = {}
    if datasets is not None:
        for dataset in datasets:
            m = re.search('^(\d+)-', os.path.basename(dataset))
            name, = m.groups()
            data = h5_utils.load(dataset)
            targets[name] = data['y'][:]

    # get average intersection percentage for each dataset
    inter = {}
    sizes = {}
    for inter_filename in inter_filenames:
        m = re.search('^(\d+)-(\d+)-', os.path.basename(inter_filename))
        a, b = m.groups()
        if a == b:
            continue  # don't count self-intersection
        data = read_pickle(inter_filename)
        if a not in inter:
            inter[a] = []
        if actives:
            assert data['inter'].shape == targets[a].shape
            sel = np.where(targets[a])[0]
            inter[a].append(np.count_nonzero(data['inter'][sel]))
            sizes[a] = np.count_nonzero(targets[a])
        else:
            inter[a].append(np.count_nonzero(data['inter']))
            sizes[a] = data['inter'].size

    for key in inter:
        inter[key] = np.mean(inter[key]) / sizes[key]

    # get scores
    df = pd.read_table(scores_filename)
    scores = {}
    for name, score in zip(df.columns[1:], df.values[-1][1:]):
        if name.startswith('PCBA'):
            name = name.split('PCBA-AID')[-1]
        elif name.startswith('MUV'):
            name = name.split('MUV-')[-1]
        elif name.startswith('TOX'):
            name = name.split('TOX-')[-1]
        elif name.startswith('DUDE'):
            name = name.split('DUDE-')[-1]
        scores[name] = score

    assert len(inter) == len(scores)
    assert np.intersect1d(inter.keys(), scores.keys()).size == len(inter)

    # plot
    x, y = [], []
    for key in inter.keys():
        x.append(inter[key])
        y.append(scores[key])
    fig = pp.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y)
    fig.savefig(output_filename, dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    args = get_args()
    main(args.inter, args.scores, args.output, args.actives, args.datasets)
