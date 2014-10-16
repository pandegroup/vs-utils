#!/usr/bin/env python
"""
Parse assay data from the NCI60 dataset.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import numpy as np

from pande_gas.utils import write_pickle
from pande_gas.utils.target_utils import Nci60Parser


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
    return parser.parse_args(input_args)


def main(input_filename, map_filename, output_filename):
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
    """
    parser = Nci60Parser(input_filename, map_filename)
    smiles, targets = parser.get_targets()

    # print the fraction of valid assay records that were found in the map
    total = np.count_nonzero(~np.isnan(
        parser.read_data(input_filename).NSC))
    print '{}/{} records matched'.format(len(targets), total)

    # save SMILES and targets
    write_pickle({'smiles': smiles, 'targets': targets}, output_filename)

if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.map, args.output)
