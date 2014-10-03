#!/usr/bin/env python
"""
Get Tox21 challenge targets from SD fields.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "3-clause BSD"

import argparse
import cPickle
import numpy as np
from rdkit_utils import serial


def get_args():
    """
    Get command-line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        help='Input filename.')
    parser.add_argument('output',
                        help='Output filename.')
    return parser.parse_args()


def main(input_filename, output_filename):
    """
    Get Tox21 chellenge targets from SD fields.

    Parameters
    ----------
    input_filename : str
        Input filename.
    output_filename : str
        Output filename.
    """
    reader = serial.MolReader()
    reader.open(input_filename)
    prop_names = ['NR-AR', 'NR-AhR', 'NR-AR-LBD', 'NR-ER', 'NR-ER-LBD',
                  'NR-Aromatase', 'NR-PPAR-gamma', 'SR-ARE', 'SR-ATAD5',
                  'SR-HSE', 'SR-MMP', 'SR-p53']
    names = []
    targets = []
    skipped = []
    for mol in reader.get_mols():
        names.append(mol.GetProp('_Name'))
        this_targets = np.ma.masked_all(len(prop_names), dtype=int)
        for prop in list(mol.GetPropNames()):
            try:
                idx = prop_names.index(prop)
                this_targets[idx] = mol.GetProp(prop)
            except ValueError:  # make sure we don't miss anything important
                if prop not in skipped:
                    skipped.append(prop)
                    print "Skipping property '{}'".format(prop)
                continue
        targets.append(this_targets)
    names = np.asarray(names)
    targets = np.ma.vstack(targets)
    assert np.all(np.sum(~targets.mask, axis=1))  # all should have a target
    with open(output_filename, 'wb') as f:
        cPickle.dump({'names': names, 'targets': targets}, f,
                     cPickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
    args = get_args()
    main(args.input, args.output)
