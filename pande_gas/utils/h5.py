"""
HDF5 file manipulation.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import h5py

save_options = {'chunks': True,
                'fletcher32': True,
                'shuffle': True,
                'compression': 'gzip',
                'compression_opts': 1}


def dump(data, filename, options=None):
    if options is None:
        options = save_options
    with h5py.File(filename) as f:
        for key, value in data.items():
            f.create_dataset(key, data=value, **options)
