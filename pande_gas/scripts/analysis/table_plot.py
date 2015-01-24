"""
Plot the data presented in Table 1.
"""
import argparse
from matplotlib import artist
from matplotlib import pyplot as pp
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import pandas as pd
import seaborn as sns
sns.set(style='whitegrid')


def get_args(input_args=None):
    """
    Get command-line arguments.

    Parameters
    ----------
    input_args : list, optional
        Input arguments. If not provided, defaults to sys.argv[1:].
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=1,
                        help='Input filename.')
    parser.add_argument('-o', '--output', required=1,
                        help='Output filename.')
    return parser.parse_args(input_args)


def main(input_filename, output_filename):
    df = pd.read_table(input_filename)
    idx = {'MUV': [], 'TOX': [], 'PCBA': [], 'DUDE': []}
    for i, c in enumerate(df.columns):
        for key in idx.iterkeys():
            if c.startswith(key):
                idx[key].append(i)
                break
    for key in idx.iterkeys():
        idx[key] = np.asarray(idx[key], dtype=int)

    data = {}
    names = []
    for model in df.values:
        name = model[0]
        names.append(name)
        for key, i in idx.iteritems():
            if key not in data:
                data[key] = []
            data[key].append(model[i])
    for key in data.iterkeys():
        data[key] = np.asarray(data[key], dtype=float).T

    fig = pp.figure(figsize=(4, 8))
    grid = ImageGrid(fig, 111, nrows_ncols=(3, 1), label_mode='L', aspect=0,
                     share_all=1, axes_pad=0.3)
    colors = sns.color_palette('hls', 8)
    pos = np.arange(len(names))
    for key, ax in zip(['PCBA', 'MUV', 'TOX'], grid):
        for i, this in enumerate(data[key].T):
            bp = ax.boxplot(this, notch=1, positions=[pos[i]], sym='o',
                            patch_artist=1, widths=0.9)
            artist.setp(bp['boxes'], facecolor=colors[i])
            artist.setp(bp['caps'], color='k')
            artist.setp(bp['whiskers'], color='k', linestyle='solid')
            artist.setp(bp['medians'], color='k')
            artist.setp(bp['fliers'], color='k', marker='o', markersize=2)
        if key == 'TOX':
            ax.set_title('Tox21')
        else:
            ax.set_title(key)
        ax.set_xticks(pos)
        ax.set_xticklabels(names, rotation=45, ha='right')
        if key == 'MUV':
            ax.set_ylabel('Median 5-Fold Average AUC')
        ax.set_xlim(-0.55, 7.55)
    fig.savefig(output_filename, dpi=300, bbox_inches='tight',
                transparent=True)

if __name__ == '__main__':
    args = get_args()
    main(args.input, args.output)
