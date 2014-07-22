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


def dump(data, filename, attrs=None, options=None):
    """
    Dump data to an HDF5 file.

    Parameters
    ----------
    data : dict
        Datasets to save.
    filename : str
        Output filename.
    attrs : dict, optional
        HDF5 attributes to set for the file.
    options : dict, optional
        Keyword arguments to create_dataset.
    """
    if options is None:
        options = save_options
    with h5py.File(filename) as f:
        for key, value in data.items():
            f.create_dataset(key, data=value, **options)
        if attrs is not None:
            for key, value in attrs.items():
                if value is None:
                    value = 'None'
                f.attrs[key] = value
