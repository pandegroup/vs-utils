"""
Plot target metrics.
"""
import argparse
import matplotlib
matplotlib.use('Agg')
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


def get_classes(filenames, skip='dude'):
    classes = {}
    corrections = {'nuclear receptor': 'transcription factor'}
    for filename in filenames:
        df = pd.read_table(filename)
        for aid, klass, subklass in zip(df['Dataset'], df['Target Class'],
                                        df['Target Subclass']):
            if aid.endswith('*'):
              aid = aid[:-1]
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

            if skip is not None and aid.startswith(skip):
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


def match(classes, scores):
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
    return data, x_labels


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
    scores, datasets = get_scores(scores_filename, False, True)  # name -> score
    data, x_labels = match(classes, scores)

    # sort by mass
    masses = np.asarray([len(a) for a in data], dtype=int)
    assert np.sum(masses) == 157, np.sum(masses)
    #means = np.asarray([np.mean(a) for a in data], dtype=float)
    #masses = means
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
    ax.set_xlabel('Target Class', fontsize=14)
    ax.set_xticklabels(x_labels[sort], rotation=30, ha='right', fontsize=14)
    #ax.set_ylabel(r'$\Delta$ Mean AUC')
    ax.set_ylabel(r'$\Delta$ Log-Odds Mean AUC', fontsize=14)
    ax.plot(x, y, '.', color='k')
    ax.set_ylim(-1, 2.5)
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

    '''
    # pie chart showing target distribution
    masses = np.asarray([len(a) for a in data], dtype=int)

    def my_autopct(pct):
        """
        See http://stackoverflow.com/questions/6170246.
        """
        total = np.sum(masses)
        val = int(pct * total / 100.)
        return '{p:.2f}% ({v:d})'.format(p=pct, v=val)

    sort = np.argsort(masses)[::-1]
    fig = pp.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax.pie(masses[sort], labels=x_labels[sort], autopct=my_autopct)
    fig.savefig('target_pie.png', dpi=300, bbox_inches='tight',
                transparent=True)
    '''
    classes = get_classes(classes_filenames, skip=None)  # class -> name
    data, x_labels = match(classes, scores)
    assert len(classes) == len(x_labels)
    masses = np.asarray([len(a) for a in data], dtype=int)
    sort = np.argsort(masses)[::-1]
    grouped_classes = np.zeros((4, len(classes)), dtype=int)
    index = ['pcba', 'muv', 'tox', 'dude']
    for j, klass in enumerate(x_labels[sort]):
        for dataset in classes[klass]:
            start = dataset.split('-')[0]
            i = index.index(start)
            grouped_classes[i, j] += 1
    print grouped_classes

    # stack everything on top of the total plot
    fig = pp.figure()
    ax = fig.add_subplot(111)
    x = np.arange(len(masses))
    colors = sns.color_palette(n_colors=len(index))
    legend_keys = ['PCBA', 'MUV', 'Tox21', 'DUD-E'][::-1]
    for i in xrange(len(grouped_classes)):
        ax.bar(x, np.sum(grouped_classes[:4-i], axis=0), color=colors[i],
               label=legend_keys[i])
    ax.set_ylabel('Count')
    ax.set_xlabel('Target Class')
    ax.set_xticks(x+0.4)
    ax.set_xticklabels(x_labels[sort], rotation=90)
    ax.grid(axis='x')
    pp.legend(loc=0)
    fig.savefig('target_bar.png', dpi=300, bbox_inches='tight',
                transparent=True)

if __name__ == '__main__':
    args = get_args()
    main(args.input, args.scores, args.output)
