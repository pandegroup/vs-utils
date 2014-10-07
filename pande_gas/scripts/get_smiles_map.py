#!/usr/bin/env python
"""
Get SMILES for compounds and map to compound names.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse

from rdkit_utils import serial

from pande_gas.utils import read_pickle, SmilesMap, write_pickle


def parse_args(input_args=None):
    """
    Parse command-line arguments.

    Parameters
    ----------
    input_args : list, optional
        Input arguments. If not provided, defaults to sys.argv[1:].
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=1, nargs='+',
                        help='Input molecule filename(s).')
    parser.add_argument('-o', '--output', required=1,
                        help='Output filename.')
    parser.add_argument('-p', '--prefix',
                        help='Prefix to prepend to molecule IDs.')
    parser.add_argument('--no-duplicates', action='store_false',
                        dest='allow_duplicates',
                        help='Whether to allow duplicate SMILES.')
    parser.add_argument('-u', '--update', action='store_true',
                        help='Update existing map with same output filename.')
    return parser.parse_args(input_args)


def main(input_filenames, output_filename, id_prefix=None,
         allow_duplicates=True, update=False):
    """
    Get SMILES for compounds and map to compound names.

    Parameters
    ----------
    input_filenames : list
        Input molecule filenames.
    output_filename : str
        Output filename.
    id_prefix : str, optional
        Prefix to prepend to IDs.
    allow_duplicates : bool, optional (default True)
        Whether to allow duplicate SMILES.
    update : bool, optional (default False)
        Whether to update an existing map with the same output filename. If
        False, a new map will be generated using only the input file(s).
    """
    smiles = SmilesMap(prefix=id_prefix, allow_duplicates=allow_duplicates)

    # update existing map
    if update:
        smiles.map = read_pickle(output_filename)

    for input_filename in input_filenames:
        print input_filename
        with serial.MolReader().open(input_filename) as reader:
            for mol in reader:
                smiles.add_mol(mol)
    write_pickle(smiles.get_map(), output_filename)

if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.output, args.prefix, args.allow_duplicates,
         args.update)
