#!/usr/bin/env python
"""
Get SMILES for compounds and map to compound names.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import cPickle
import gzip

from rdkit_utils import serial

from pande_gas.utils import SmilesMap


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
    parser.add_argument('-u', '--update', action='store_true',
                        help='Update existing map with same output filename.')
    return parser.parse_args(input_args)


def main(input_filenames, output_filename, id_prefix=None, update=False):
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
    update : bool, optional (default False)
        Whether to update an existing map with the same output filename. If
        False, a new map will be generated using only the input file(s).
    """
    smiles = SmilesMap(id_prefix)

    # update existing map
    if update:
        if output_filename.endswith('.gz'):
            f = gzip.open(output_filename)
        else:
            f = open(output_filename)
        current_map = cPickle.load(f)
        f.close()
        smiles.map = current_map

    for input_filename in input_filenames:
        print input_filename
        with serial.MolReader().open(input_filename) as reader:
            for mol in reader:
                smiles.add_mol(mol)
    if output_filename.endswith('.gz'):
        f = gzip.open(output_filename, 'wb')
    else:
        f = open(output_filename, 'wb')
    cPickle.dump(smiles.get_map(), f, cPickle.HIGHEST_PROTOCOL)
    f.close()

if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.output, args.prefix, args.update)
