"""
Utilities for working with sklearn.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import numpy as np

from sklearn.cross_validation import LeavePLabelOut
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import LinearSVC, SVC, SVR


class RandomLeavePLabelOut(object):
    """
    Randomly choose a subset of possible LeavePLabelOut experiments (without
    replacement).

    Parameters
    ----------
    labels : array_like
        Labels used for partitioning.
    p : int
        Number of labels to leave out of training data in each split.
    n_iter : int, optional (default 10)
        Number of experiments to perform.
    random_state : int, optional
        Random state.
    """
    def __init__(self, labels, p, n_iter=10, random_state=None):
        self.cv = LeavePLabelOut(labels, int(p))
        self.rng = np.random.RandomState(random_state)
        self.splits = self.rng.shuffle(list(self.cv))
        self.n_iter = n_iter
        self.random_state = random_state

    def __iter__(self):
        for i, (train, test) in enumerate(self.splits):
            if i >= self.n_iter:
                break
            yield train, test


class DefaultGridSearch(object):
    """
    Default grid search parameters for several estimator types.

    Reference
    ---------
    Chih-Wei Hsu, Chih-Chung Chang, and Chih-Jen Lin,
    A Practical Guide to Support Vector Classification.
    http://www.csie.ntu.edu.tw/~cjlin/papers/guide/guide.pdf

    Parameters
    ----------
    estimator : object
        sklearn estimator.
    """
    def __init__(self, estimator):
        self.estimator = estimator

    def get_params(self):
        """
        Get default grid search parameters.
        """
        params = False
        if isinstance(self.estimator, GaussianNB):
            params = None
        elif isinstance(self.estimator, LinearSVC):
            params = {'C': self.get_default_C_params()}
        elif isinstance(self.estimator, LogisticRegression):
            params = {'C': self.get_default_C_params()}
        elif isinstance(self.estimator, RandomForestClassifier):
            params = {'criterion': ['gini', 'entropy'],
                      'max_features': ['sqrt', 'log2']}
        elif isinstance(self.estimator, SVC):
            if self.estimator.kernel == 'linear':
                params = {'C': self.get_default_C_params()}
            elif self.estimator.kernel == 'rbf':
                params = {'C': self.get_default_C_params(),
                          'gamma': self.get_default_gamma_params()}
        elif isinstance(self.estimator, SVR):
            if self.estimator.kernel == 'rbf':
                params = {'C': self.get_default_C_params(),
                          'epsilon': self.get_default_epsilon_params(),
                          'gamma': self.get_default_gamma_params()}
        if params is False:
            raise NotImplementedError('Estimator {} is not supported.'.format(
                self.estimator))
        return params

    @staticmethod
    def get_default_C_params():
        """
        Get default parameters for C.
        """
        return np.logspace(-5, 15, 11, base=2.0)

    @staticmethod
    def get_default_gamma_params():
        """
        Get default parameters for gamma.
        """
        return np.logspace(-15, 3, 10, base=2.0)

    @staticmethod
    def get_default_epsilon_params():
        """
        Get default parameters for epsilon.
        """
        return np.logspace(-9, 0, 10, base=2.0)


def resolve_estimator(name):
    """
    Resolve the estimator type from a string.

    Parameters
    ----------
    name : str
        Estimator name.
    """
    if name == 'lr':
        estimator = LogisticRegression()
    elif name == 'svm':
        estimator = SVC()
    elif name == 'lsvm':
        estimator = LinearSVC()
    elif name == 'nb':
        estimator = GaussianNB()
    elif name == 'rf':
        estimator = RandomForestClassifier(1000)
    else:
        raise NotImplementedError(name)
    return estimator
