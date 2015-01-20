"""
Plot target metrics.
"""
import argparse
import matplotlib.pyplot as pp
import numpy as np
import pandas as pd
from scipy.stats import linregress
import seaborn as sns
sns.set(style='whitegrid')

from pande_gas.scripts.analysis import get_scores


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
        for aid, klass, subklass in zip(df['Dataset'], df['Target Class'],
                                        df['Target Subclass']):
            if 'Bad' in df.columns:
                marks = df['Bad'][df['Dataset'] == aid].values
                assert len(marks) == 1
                if marks[0] == 'x':
                    print 'SKIPPING', aid
                    continue

            # handle subclasses
            if klass.lower() == 'enzyme':
                if (isinstance(subklass, str) and
                        subklass.lower() in ['protein kinase', 'protease']):
                    klass = subklass
                else:
                    klass = 'other enzyme'

            # class corrections
            if klass in corrections:
                klass = corrections[klass]

            if klass not in classes:
                classes[klass] = []

            if aid.startswith('dude'):
                print 'SKIPPING', aid
                continue  # skipping dude

            classes[klass].append(str(aid).lower())

    # merge small classes into misc.
    merge = []
    for klass, members in classes.iteritems():
        if len(members) < 5:
            merge.append(klass)
    for klass in merge:
        print 'Merging {} into misc.'.format(klass)
        for member in classes[klass]:
            classes['miscellaneous'].append(member)
        del classes[klass]

    return classes


def main(classes_filenames, scores_filename, output_filename):
    """

    Notes:
        * Count 'nuclear receptor' as 'transcription factor'.
        * Subdivide enzymes into
            - protein kinases
            - proteases
            - other enzymes
    """
    classes = get_classes(classes_filenames)  # class -> name
    scores, datasets = get_scores(scores_filename, False)  # name -> score
    data = []
    x_labels = []
    for klass in classes.keys():
        x_labels.append(klass)
        class_data = []
        for name in classes[klass]:
            if name.startswith('tox'):
                name += '-train'
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
    assert np.sum(masses) == 157, np.sum(masses)
    means = np.asarray([np.mean(a) for a in data], dtype=float)
    masses = means
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
    fig.savefig(output_filename, dpi=300, bbox_inches='tight',
                transparent=True)

    # correlation
    x, y = [], []
    for mass, this in zip(masses[sort], data[sort]):
        for score in this:
            x.append(mass - 1)  # count number of other datasets in this class
            y.append(score)
    m, b, r, p, s = linregress(x, y)
    print r, r ** 2

    # pie chart showing target distribution
    masses = np.asarray([len(a) for a in data], dtype=int)
    sort = np.argsort(masses)[::-1]
    fig = pp.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax.pie(masses[sort], labels=x_labels[sort], autopct='%d')
    fig.savefig('target_pie.png', dpi=300, bbox_inches='tight',
                transparent=True)

if __name__ == '__main__':
    args = get_args()
    main(args.input, args.scores, args.output)
