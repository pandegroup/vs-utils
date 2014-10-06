#!/usr/bin/env python
"""
Get Tox21 challenge datasets.

Notes:
* Some SMILES strings are represented by multiple compounds, with some overlap
of assays. These compounds need to be condensed and their assay outcomes need
to be reconciled.
* For compounds with activities that do not agree...look at the assays if
possible.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "3-clause BSD"

import argparse
import cPickle
import gzip
import numpy as np

from rdkit import Chem

from rdkit_utils import serial


def get_args():
    """
    Get command-line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        help='Input filename.')
    parser.add_argument('--merge', choices=['max', 'min', 'majority-',
                                            'majority+'],
                        help='Target merge strategy.')
    return parser.parse_args()


def main(input_filename, merge):
    """
    Get Tox21 chellenge datasets.

    Parameters
    ----------
    input_filename : str
        Input filename.
    merge : str
        Target merge strategy.
    """
    reader = serial.MolReader()
    reader.open(input_filename)
    dataset_names = ['NR-AR', 'NR-AhR', 'NR-AR-LBD', 'NR-ER', 'NR-ER-LBD',
                     'NR-Aromatase', 'NR-PPAR-gamma', 'SR-ARE', 'SR-ATAD5',
                     'SR-HSE', 'SR-MMP', 'SR-p53']
    data = {dataset: {} for dataset in dataset_names}
    skipped = []
    for mol in reader.get_mols():
        smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True,
                                  canonical=True)
        for prop in list(mol.GetPropNames()):
            if prop in dataset_names:
                score = int(mol.GetProp(prop))
                if smiles not in data[prop]:
                    data[prop][smiles] = []
                data[prop][smiles].append(score)
            else:  # make sure we don't miss anything important
                if prop not in skipped:
                    skipped.append(prop)
                    print "Skipping property '{}'".format(prop)
                continue

    # reconcile activity differences
    for dataset in dataset_names:
        for smiles, targets in data[dataset].items():
            targets = np.asarray(targets, dtype=int)
            if len(targets) == 1:
                data[dataset][smiles] = targets[0]
            elif np.all(targets == 0):
                data[dataset][smiles] = 0
            elif np.all(targets == 1):
                data[dataset][smiles] = 1
            else:
                if merge == 'max':
                    data[dataset][smiles] = max(targets)
                elif merge == 'min':
                    data[dataset][smiles] = min(targets)
                elif merge == 'majority-':  # 0.5 rounds down
                    data[dataset][smiles] = int(np.round(np.mean(targets)))
                elif merge == 'majority+':  # 0.5 rounds up
                    data[dataset][smiles] = (int(np.round(
                        np.mean(targets) + 1)) - 1)

    # save individual datasets
    for dataset in dataset_names:
        assert len(data[dataset])  # check for at least one molecule

        # build mol and target arrays
        mols = []
        targets = []
        for mol, target in data[dataset].items():
            mols.append(mol)
            targets.append(target)
        targets = np.asarray(targets, dtype=int)
        pos = np.count_nonzero(targets == 1)
        neg = np.count_nonzero(targets == 0)
        assert pos + neg == targets.size
        print '{}\t{}\t{}'.format(dataset, pos, neg)
        with gzip.open('{}.ism.gz'.format(dataset), 'wb') as f:
            for i, mol in enumerate(mols):
                f.write('{}\t{}\n'.format(mol, '{}-{}'.format(dataset, i)))
        with open('{}-targets.pkl'.format(dataset), 'wb') as f:
            cPickle.dump(targets, f, cPickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
    args = get_args()
    main(args.input, args.merge)
