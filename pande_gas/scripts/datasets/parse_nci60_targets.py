#!/usr/bin/env python
"""
Parse assay data from the NCI60 dataset and write a separate target file for
each dataset.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import numpy as np
import os

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
    parser.add_argument('-d', '--dir', default='.',
                        help='Directory in which to write target files.')
    parser.add_argument('-p', '--prefix', default='nci60',
                        help='Prefix for target files.')
    parser.add_argument('-s', '--suffix', default='pkl.gz',
                        choices=['pkl', 'pkl.gz'],
                        help='Suffix for target files.')
    return parser.parse_args(input_args)


def main(input_filename, map_filename, directory='.', prefix='nci60',
         suffix='pkl.gz'):
    """
    Get regression targets.

    Parameters
    ----------
    input_filename : str
        PCBA data filename.
    map_filename : str
        ID->SMILES map filename.
    directory : str, optional (default '.')
        Directory in which to write target files.
    prefix : str, optional (default 'nci60')
        Prefix for target files.
    suffix : str, optional (default 'pkl.gz')
        Suffix for target files.
    """
    parser = Nci60Parser(input_filename, map_filename)
    split_targets = parser.split_targets()

    # get total record count
    total = np.count_nonzero(~np.isnan(parser.read_data().NSC))

    # write a separate file for each dataset
    # note that split_targets is an OrderedDict
    for i, name in enumerate(split_targets.keys()):
        data = split_targets[name]
        # print the fraction of valid assay records that were found in the map
        print '{}\t{}/{} records matched'.format(
            name, len(data['targets']), total)
        write_pickle(
            data,
            os.path.join(directory,
                         '{}-{:02}-targets.{}'.format(prefix, i, suffix)))

if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.map, args.dir, args.prefix, args.suffix)
