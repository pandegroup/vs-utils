"""
Plot target metrics.
"""
import argparse
import matplotlib.pyplot as pp
import numpy as np
import pandas as pd
from scipy.stats import linregress
import seaborn as sns

from pande_gas.utils import read_pickle


def get_args(input_args=None):
    """
    Get command-line arguments.

    Parameters
    ----------
    input_args : list, optional
        Input arguments. If not provided, defaults to sys.argv[1:].
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=1, nargs='+',
                        help='Target classes.')
    parser.add_argument('-s', '--scores', required=1,
                        help='Scores.')
    parser.add_argument('-o', '--output', required=1,
                        help='Output filename.')
    return parser.parse_args(input_args)


def get_classes(filenames):
    classes = {}
    corrections = {'nuclear receptor': 'transcription factor'}
    for filename in filenames:
        df = pd.read_table(filename)
        for aid, klass in zip(df['AID'], df['Target Class']):
            if 'Bad' in df.columns:
                marks = df['Bad'][df['AID'] == aid].values
                assert len(marks) == 1
                if marks[0] == 'x':
                    print 'SKIPPING', aid
                    continue

            # class corrections
            if klass in corrections:
                klass = corrections[klass]

            if klass not in classes:
                classes[klass] = []
            classes[klass].append(str(aid))
    return classes


def get_scores(filename):
    df = pd.read_table(filename)
    scores = {}
    ref_idx = 1
    new_idx = 27
    print df.values[ref_idx][0], df.values[new_idx][0]  # print scores
    datasets = {'PCBA': [], 'MUV': [], 'TOX': [], 'DUDE': []}
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
    return scores, datasets


def main(classes_filenames, scores_filename, output_filename):
    """

    Notes:
        * Count 'nuclear receptor' as 'transcription factor'.
    """
    classes = get_classes(classes_filenames)  # class -> name
    scores, datasets = get_scores(scores_filename)  # name -> score
    data = []
    x_labels = []
    for klass in classes.keys():
        x_labels.append(klass)
        class_data = []
        for name in classes[klass]:
            try:
                class_data.append(scores[name])
            except KeyError as e:
                import IPython
                IPython.embed()
                raise e
        data.append(class_data)
    data = np.asarray(data)
    x_labels = np.asarray(x_labels)

    # sort by mass
    masses = np.asarray([len(a) for a in data], dtype=int)
    sort = np.argsort(masses)[::-1]
    for klass, mass in zip(x_labels[sort], masses[sort]):
        print klass, mass

    # get x and y for manual points
    x, y = [], []
    for i, this in enumerate(data[sort]):
        for score in this:
            x.append(i + 1)
            y.append(score)

    # plot
    fig = pp.figure()
    ax = fig.add_subplot(111)
    #sns.violinplot(data[sort], inner='points', ax=ax, positions=masses[sort])
    #ax.set_xticks(masses[sort])
    sns.violinplot(data[sort], inner='points', ax=ax)
    ax.set_xlabel('Target Class')
    ax.set_xticklabels(x_labels[sort], rotation=90)
    ax.set_ylabel(r'$\Delta$ Mean AUC')
    ax.plot(x, y, '.', color='k')
    fig.savefig(output_filename, dpi=300, bbox_inches='tight')

    # correlation
    x, y = [], []
    for mass, this in zip(masses[sort], data[sort]):
        for score in this:
            x.append(mass - 1)  # count number of other datasets in this class
            y.append(score)
    m, b, r, p, s = linregress(x, y)
    print r, r ** 2
    fig = pp.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y, '.', color='k')
    ax.set_xlabel('Number of Like Targets')
    ax.set_xticks(masses[sort] - 1)
    #ax.set_xticklabels(x_labels[sort], rotation=90)
    ax.set_ylabel(r'$\Delta$ Mean AUC')
    ax.plot([0, 136], [b, 137*m + b])
    fig.savefig('scaled.png', dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    args = get_args()
    main(args.input, args.scores, args.output)
