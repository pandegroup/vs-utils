#!/usr/bin/env python
"""
Split a dataset into chunks.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse

from vs_utils.utils import DatasetSharder


def parse_args(input_args=None):
    """
    Parse command-line arguments.

    Parameters
    ----------
    input_args : list, optional
        Input arguments. If not provided, defaults to sys.argv[1:].
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        help='Input filename.')
    parser.add_argument('-n', type=int, default=1000,
                        help='Number of molecules per chunk.')
    parser.add_argument('-p', '--prefix',
                        help='Prefix for output files. Defaults to prefix ' +
                             'of input filename.')
    parser.add_argument('-f', '--flavor', default='pkl.gz',
                        help='Output flavor.')
    return parser.parse_args(input_args)


def main(filename, shard_size, prefix, flavor):
    """
    Split a dataset into chunks.

    Parameters
    ----------
    filename : str
        Input filename.
    shard_size : int
        Number of molecules per shard.
    prefix : str
        Prefix for output files. Defaults to the prefix of the input filename,
        as extracted by guess_prefix.
    flavor : str
        Output molecule format used as the extension for shard filenames.
    """
    sharder = DatasetSharder(filename=filename, shard_size=shard_size,
                             prefix=prefix, flavor=flavor)
    sharder.shard()


if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.n, args.prefix, args.flavor)
