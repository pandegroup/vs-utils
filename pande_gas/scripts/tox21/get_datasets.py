#!/usr/bin/env python
"""
Get Tox21 challenge datasets.
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
    return parser.parse_args()


def main(input_filename):
    """
    Get Tox21 chellenge datasets.

    Parameters
    ----------
    input_filename : str
        Input filename.
    """
    reader = serial.MolReader()
    reader.open(input_filename)
    dataset_names = ['NR-AR', 'NR-AhR', 'NR-AR-LBD', 'NR-ER', 'NR-ER-LBD',
                     'NR-Aromatase', 'NR-PPAR-gamma', 'SR-ARE', 'SR-ATAD5',
                     'SR-HSE', 'SR-MMP', 'SR-p53']
    mols = {dataset: [] for dataset in dataset_names}
    targets = {dataset: [] for dataset in dataset_names}
    skipped = []
    for mol in reader.get_mols():
        for prop in list(mol.GetPropNames()):
            if prop in dataset_names:
                mols[prop].append(mol)
                targets[prop].append(int(mol.GetProp(prop)))
            else:  # make sure we don't miss anything important
                if prop not in skipped:
                    skipped.append(prop)
                    print "Skipping property '{}'".format(prop)
                continue

    # check that each dataset has at least one molecule
    for dataset in dataset_names:
        assert len(mols[dataset])
        assert len(mols[dataset]) == len(targets[dataset])
        targets[dataset] = np.asarray(targets[dataset], dtype=int)
        print '{}\t{}\t{}'.format(dataset,
                                  np.count_nonzero(targets[dataset] == 1),
                                  np.count_nonzero(targets[dataset] == 0))

    # save individual datasets
    for dataset in dataset_names:
        writer = serial.MolWriter()
        writer.open('{}.sdf.gz'.format(dataset))
        writer.write(mols[dataset])
        with open('{}-targets.pkl'.format(dataset), 'wb') as f:
            cPickle.dump(targets[dataset], f, cPickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
    args = get_args()
    main(args.input)
