"""
Tests for miscellaneous utilities.
"""
import numpy as np
import unittest

from pande_gas.utils import pad_array


class TestUtils(unittest.TestCase):
    """
    Tests for miscellaneous utilities.
    """
    def test_pad_matrix(self):
        """Pad matrix."""
        x = np.random.random((5, 6))
        assert pad_array(x, (10, 12)).shape == (10, 12)
        assert pad_array(x, 10).shape == (10, 10)
