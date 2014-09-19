#!/usr/bin/env python
"""
Build an sklearn model using pande-gas features.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "3-clause BSD"

import argparse
import cPickle
import gzip
import numpy as np
from scipy import optimize

from sklearn.grid_search import GridSearchCV
from sklearn.metrics import roc_auc_score

from pande_gas.utils import sklearn_utils


def get_args():
    """
    Get command-line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--train', required=1,
                        help='Filename containing training set features and ' +
                             'targets.')
    parser.add_argument('--test', required=1,
                        help='Filename containing test set features.')
    parser.add_argument('-e', '--estimator', required=1)
    parser.add_argument('-n', '--n_iter', type=int, default=10)
    parser.add_argument('--test-size', type=float, default=0.1)
    parser.add_argument('--seed', type=int)
    parser.add_argument('-p', '--prefix', required=1)
    parser.add_argument('-np', '--n-jobs', type=int, default=1)
    return parser.parse_args()


def read_training_data(filename):
    """
    Read training data.

    Parameters
    ----------
    filename : str
        Training set (features and labels) filename.
    """
    if filename.endswith('.gz'):
        f = gzip.open(filename)
    else:
        f = open(filename)
    data = cPickle.load(f)
    f.close()
    return data['features'], data['y'], data['scaffolds']


def read_test_data(filename):
    """
    Read test data.

    Parameters
    ----------
    filename : str
        Test set (features) filename.
    """
    if filename.endswith('.gz'):
        f = gzip.open(filename)
    else:
        f = open(filename)
    data = cPickle.load(f)
    f.close()
    return data['features'], data['names']


def save_predictions(y_pred, names, threshold, filename):
    """
    Save predictions.

    Parameters
    ----------
    y_pred : array_like
        Predicted positive class probabilities.
    names : array_like
        Compound names.
    threshold : float
        Probability threshold for classification.
    filename : str
        Output filename.
    """
    labels = np.asarray(y_pred > threshold, dtype=int)
    with open(filename, 'wb') as f:
        f.write('Sample ID\tScore\tActivity\n')
        for name, pred, label in zip(names, y_pred, labels):
            f.write('{}\t{}\t{}\n'.format(name, pred, label))


def save_model(estimator, filename):
    """
    Save model.

    Parameters
    ----------
    estimator : object
        sklearn estimator.
    filename : str
        Output filename.
    """
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'wb')
    else:
        f = open(filename, 'wb')
    cPickle.dump(estimator, f, cPickle.HIGHEST_PROTOCOL)
    f.close()


def guess_threshold(estimator, X, y):
    """
    Choose classification threshold that gives best AUC.
    """
    x0 = [np.mean(estimator.predict(X))]
    rval = optimize.minimize(_guess_threshold, x0, (estimator, X, y),
                             options={'disp': True})
    return rval.x


def _guess_threshold(threshold, estimator, X, y):
    """
    Calculate negative ROC-AUC using this classification threshold.
    """
    rval = roc_auc_score(
        y, np.asarray(estimator.predict(X) > threshold, dtype=int))
    return -rval


def main(train_filename, test_filename, estimator, output_prefix, n_iter=10,
         test_size=0.1, random_state=None, n_jobs=1):
    """
    Fit an sklearn model and predict test set labels.

    Parameters
    ----------
    train_filename : str
        Training set (features and labels) filename.
    test_filename : str
        Test set (features) filename.
    estimator : str
        Estimator name.
    output_prefix : str
        Prefix for output files.
    n_iter : int, optional (default 10)
        Number of cross-validation experiments.
    test_size : float, optional (default 0.1)
        Fraction of scaffolds in cross-validation test sets.
    random_state : int or RandomState, optional
        Random state.
    n_jobs : int, optional (default 1)
        Number of jobs to use for GridSearchCV.
    """
    estimator = sklearn_utils.resolve_estimator(estimator)
    X_train, y_train, labels = read_training_data(train_filename)
    X_test, names_test = read_test_data(test_filename)
    params = sklearn_utils.DefaultGridSearch(estimator).get_params()
    p = np.ceil(np.unique(labels).size * test_size)
    cv = sklearn_utils.RandomLeavePLabelOut(labels, p, n_iter, random_state)
    engine = GridSearchCV(estimator, params, cv=cv, scoring='roc-auc',
                          n_jobs=n_jobs)
    engine.fit(X_train, y_train)
    threshold = guess_threshold(engine, X_train, y_train)
    y_pred = engine.predict(X_test)
    save_predictions(y_pred, names_test, threshold,
                     '{}-pred.txt'.format(output_prefix))
    save_model(engine.best_estimator_, '{}-model.pkl.gz'.format(output_prefix))

if __name__ == '__main__':
    args = get_args()
    main(args.train, args.test, args.estimator, args.prefix, args.n_iter,
         args.test_size, args.seed, args.n_jobs)
