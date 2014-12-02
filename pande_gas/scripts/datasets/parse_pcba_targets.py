#!/usr/bin/env python
"""
Parse assay data from PubChem BioAssay (PCBA).
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import numpy as np

from pande_gas.utils import write_pickle
from pande_gas.utils.target_utils import PcbaParser


def parse_args(input_args=None):
    """
    Parse command-line arguments.

    Parameters
    ----------
    input_args : list, optional
        Input arguments. If not provided, defaults to sys.argv[1:].
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=1,
                        help='Input data file.')
    parser.add_argument('-m', '--map', required=1,
                        help='Molecule ID to SMILES map.')
    parser.add_argument('-o', '--output', required=1,
                        help='Output filename.')
    parser.add_argument('-c', '--cols', nargs='+', type=int,
                        help='Column indices to include.')
    return parser.parse_args(input_args)


def main(input_filename, map_filename, output_filename, column_indices=None):
    """
    Get regression targets.

    Parameters
    ----------
    input_filename : str
        PCBA data filename.
    map_filename : str
        ID->SMILES map filename.
    output_filename : str
        Output filename.
    column_indices : list, optional
        Column indices to include. If None, compounds are classified by
        activity.
    """
    parser = PcbaParser(input_filename, map_filename,
                        column_indices=column_indices)
    if column_indices is not None:
        print "Extracting data from the following columns:"
        for col in parser.get_column_names():
            print '\t', col
    smiles, targets = parser.get_targets()

    # print the fraction of valid assay records that were found in the map
    total = np.count_nonzero(~np.isnan(parser.read_data().PUBCHEM_CID))
    print '{}/{} records matched'.format(len(targets), total)
    print '\t'.join(['{}:{}'.format(i, np.count_nonzero(targets == i))
                     for i in np.unique(targets)])

    # save SMILES and targets
    write_pickle({'smiles': smiles, 'targets': targets}, output_filename)

if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.map, args.output, args.cols)
