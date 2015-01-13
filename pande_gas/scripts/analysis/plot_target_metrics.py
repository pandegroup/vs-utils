"""
Plot target metrics.
"""
import argparse
import matplotlib.pyplot as pp
import numpy as np
import pandas as pd
from scipy.stats import linregress
import seaborn as sns

from . import get_scores


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


if __name__ == '__main__':
    args = get_args()
    main(args.input, args.scores, args.output)
