"""
Tests for HDF5 utilities.
"""
import h5py
import numpy as np
import os
import tempfile

from pande_gas.utils import h5_utils


def test_dump():
    """Test dump to HDF5."""
    _, filename = tempfile.mkstemp()

    # write some data
    data = {}
    data['a'] = np.random.random((10, 3))
    data['b'] = np.random.randint(6, size=10)
    h5_utils.dump(data, filename)

    # make sure we can read it
    with h5py.File(filename) as f:
        assert np.array_equal(data['a'], f['a'])
        assert np.array_equal(data['b'], f['b'])

    # cleanup
    os.remove(filename)
