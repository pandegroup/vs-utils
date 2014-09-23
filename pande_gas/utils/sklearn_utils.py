"""
Utilities for working with sklearn.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

from collections import OrderedDict
import numpy as np

from sklearn.cross_validation import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import LinearSVC, SVC, SVR


class StratifiedClusterKFold(object):
    """
    Cluster-based partitioning with stratification according to cluster labels.

    Parameters
    ----------
    y : array_like
        Class labels for compounds.
    clusters : array_like
        Clusters for compounds.
    n_folds : int, optional (default 3)
        Number of folds.
    shuffle : bool, optional (default False)
        Whether to shuffle before stratification.
    random_state : int or RandomState, optional
        Random state used for shuffling.
    """
    def __init__(self, y, clusters, n_folds=3, shuffle=False,
                 random_state=None):
        self.y = y
        self.clusters = clusters
        self.cluster_labels = self.get_cluster_labels()
        self.cv = StratifiedKFold(self.cluster_labels.values(),
                                  n_folds=n_folds, shuffle=shuffle,
                                  random_state=random_state)

    def __iter__(self):
        """
        Yield compounds in clusters that are separated with StratifiedKFold.
        """
        for cluster_train, cluster_test in self.cv:
            train = self.get_compound_indices(cluster_train)
            test = self.get_compound_indices(cluster_test)
            yield train, test

    def get_compound_indices(self, cluster_indices):
        """
        Get compound indices matching cluster indices.

        Parameters
        ----------
        cluster_indices : array_like
            Indices of selected clusters.
        """
        this_clusters = np.asarray(self.cluster_labels.keys())[cluster_indices]
        mask = np.zeros_like(self.y, dtype=bool)
        for cluster in np.atleast_1d(this_clusters):
            mask[self.clusters == cluster] = True
        return np.where(mask)[0]

    def get_cluster_labels(self):
        """
        Assign labels to clusters based on compound labels.

        Clusters are marked as either 'active' or 'inactive' depending on
        whether they contain _any_ active compounds.
        """
        cluster_labels = OrderedDict()
        for cluster in np.unique(self.clusters):
            if np.any(self.y[self.clusters == cluster]):
                cluster_labels[cluster] = 1
            else:
                cluster_labels[cluster] = 0
        return cluster_labels

    def __len__(self):
        return len(self.cv)


class RandomLeavePLabelOut(object):
    """
    Randomly choose a subset of possible LeavePLabelOut experiments (without
    replacement). It takes too long to enumerate all possible experiments and
    select from that list, so we instead choose up to n_iter experiments.

    Parameters
    ----------
    labels : array_like
        Labels used for partitioning.
    p : int
        Number of labels to leave out of training data in each split.
    n_iter : int, optional (default 10)
        Number of experiments to perform.
    max_attempts : int, optional (default 100)
        Maximum number of attempts to sample a new experiment.
    random_state : int, optional
        Random state.
    """
    def __init__(self, labels, p, n_iter=10, max_attempts=100,
                 random_state=None):
        self.labels = labels
        self.unique_labels, self.unique_inverse = np.unique(
            labels, return_inverse=True)
        assert p > 0 and p == int(p), 'p must be a positive integer.'
        self.p = p
        self.rng = np.random.RandomState(random_state)
        self.max_attempts = max_attempts
        self.n_iter = n_iter

    def __iter__(self):
        """
        Choose up to self.n_iter experiments, failing if a new experiment
        cannot be generated after self.max_attempts attempts.
        """
        states = []
        n_attempts = 0
        while len(states) < self.n_iter and n_attempts < self.max_attempts:
            state = self.rng.get_state()[1][0]
            if state in states:
                n_attempts += 1
                continue
            else:
                n_attempts = 0
                states.append(state)
            labels = self.unique_labels.copy()
            self.rng.shuffle(labels)
            train_labels = labels[:-self.p]
            test_labels = labels[-self.p:]
            assert not np.intersect1d(train_labels, test_labels).size
            train = self._get_indices(train_labels)
            test = self._get_indices(test_labels)
            assert train.size + test.size == self.labels.size
            assert not np.intersect1d(train, test).size
            yield train, test

    def __len__(self):
        return self.n_iter

    def _get_indices(self, which_labels):
        """
        Get the indices of self.labels that match labels in which_labels.

        Parameters
        ----------
        which_labels : array_like
            Labels to select.
        """
        indices = np.zeros_like(self.labels, dtype=bool)
        for label in which_labels:
            indices[self.labels == label] = True
        return np.where(indices)[0]  # convert mask to indices


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
