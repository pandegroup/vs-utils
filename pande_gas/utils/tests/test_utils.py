"""
Tests for miscellaneous utilties.
"""
import numpy as np

from pande_gas.utils import pad_array


def test_pad_matrix():
    """Pad matrix."""
    x = np.random.random((5, 6))
    assert pad_array(x, (10, 12)).shape == (10, 12)
    assert pad_array(x, 10).shape == (10, 10)
