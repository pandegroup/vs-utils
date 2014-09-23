"""
Tests for sklearn_utils.
"""

import numpy as np
import unittest

from .. import sklearn_utils


class TestStratifiedClusterKFold(unittest.TestCase):
    """
    Test StratifiedClusterKFold.
    """
    def setUp(self):
        """
        Set up tests.
        """
        self.y = np.random.randint(2, size=100)  # two classes
        self.clusters = np.random.randint(10, size=100)  # ten clusters
        self.cv = sklearn_utils.StratifiedClusterKFold(self.y, self.clusters)

    def test_get_cluster_labels(self):
        """
        Test StratifiedClusterKFold.get_cluster_labels.

        Verify that clusters containing active compounds are marked active.
        """
        cluster_labels = self.cv.get_cluster_labels()
        for cluster, label in cluster_labels.items():
            this_y = self.y[self.clusters == cluster]
            if np.any(this_y):
                assert label
            else:
                assert not label

    def test_get_compound_indices(self):
        """
        Test StratifiedClusterKFold.get_compound_indices.

        Verify that cluster indices are correctly mapped to compound indices.
        """

        # individual clusters
        for i, cluster in enumerate(self.cv.cluster_labels.keys()):
            indices = np.where(self.clusters == cluster)[0]
            assert np.array_equal(self.cv.get_compound_indices(i), indices)

        # all clusters
        all_cluster_indices = np.arange(len(self.cv.cluster_labels))
        assert np.array_equal(
            self.cv.get_compound_indices(all_cluster_indices),
            np.arange(len(self.y)))

        # no clusters
        assert np.array_equal(self.cv.get_compound_indices([]), [])

    def verify_cv(self, cv):
        """
        Test proper fold generation.

        Parameters
        ----------
        cv : StratifiedClusterKFold
            Cross-validation subset iterator.
        """
        for train, test in cv:

            # basic tests
            assert train.size + test.size == self.y.size  # full coverage
            assert not np.intersect1d(train, test).size  # no leaks
            assert np.unique(train).size == train.size  # no duplication
            assert np.unique(test).size == test.size  # no duplication

            # cluster-based tests
            # check that each cluster is in either train or test only
            for cluster in np.unique(self.clusters):
                this_cluster = np.where(self.clusters == cluster)[0]
                all_in_train = np.all(np.in1d(this_cluster, train))
                all_in_test = np.all(np.in1d(this_cluster, test))
                assert all_in_train or all_in_test
                assert not (all_in_train and all_in_test)

    def test_default(self):
        """
        Test StratifiedClusterKFold with default parameters.
        """
        self.verify_cv(self.cv)

    def test_shuffled(self):
        """
        Test StratifiedClusterKFold with shuffling.
        """
        self.cv = sklearn_utils.StratifiedClusterKFold(self.y, self.clusters,
                                                       shuffle=True)
        self.verify_cv(self.cv)
