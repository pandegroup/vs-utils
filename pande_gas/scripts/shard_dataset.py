#!/usr/bin/env python
"""
Split a dataset into chunks.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import os

from rdkit_utils import PicklableMol, serial


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
    parser.add_argument('-f', '--format', default='pkl.gz',
                        help='Output format.')
    return parser.parse_args(input_args)


def main(input_filename, chunk_size=1000, output_prefix=None,
         output_format='pkl.gz'):
    """
    Split a dataset into chunks.

    Parameters
    ----------
    input_filename : str
        Input filename.
    chunk_size : int, optional (default 1000)
        Number of molecules per chunk.
    output_prefix : str, optional
        Prefix for output files. Defaults to the prefix of input_filename.
    output_format : str, optional (default 'pkl.gz')
        Extension for output files.
    """
    if output_prefix is None:
        output_prefix = os.path.basename(input_filename).split('.')[0]
    reader = serial.MolReader()
    reader.open(input_filename)
    filename_generator = get_filenames(output_prefix, output_format)
    mols = []
    for mol in reader.get_mols():
        mols.append(mol)
        if len(mols) >= chunk_size:
            write(mols, filename_generator.next())
            mols = []
    if len(mols):
        write(mols, filename_generator.next())


def get_filenames(prefix, extension, index=0):
    """
    Generate shard filenames.

    Parameters
    ----------
    prefix : str
        Prefix for filenames.
    extension : str
        Extension for filenames.
    index : int, optional (default 0)
        Initial shard index.
    """
    while True:
        filename = '{}-{}.{}'.format(prefix, index, extension)
        yield filename
        index += 1


def write(mols, filename):
    """
    Write molecules to file.

    Molecules are converted to PicklableMol to preserve properties, such as
    molecule names.

    Parameters
    ----------
    mols : array_like
        Molecules.
    filename : str
        Output filename.
    """
    mols = [PicklableMol(mol) for mol in mols]  # preserve properties
    writer = serial.MolWriter()
    print '{}: {}'.format(filename, len(mols))
    with writer.open(filename) as f:
        f.write(mols)


if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.n, args.prefix, args.format)
