#!/usr/bin/env python
"""
Generate conformers with the RDKit.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
from joblib import delayed, Parallel

from pande_gas.utils import rdkit_utils as rd


def main():
    confs = Parallel(n_jobs=args.n_jobs, verbose=5)(
        delayed(rd.generate_conformers)(mol, args.n_conformers,
                                        args.rmsd_threshold)
        for mol in rd.read_mols(args.input))
    passed = []
    for i, mol in enumerate(confs):
        if mol is None:
            print '{} failed.'.format(mol.GetProp('_Name'))
            continue
        passed.append(mol)
    rd.write_mols(passed, args.output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        help='Input SDF molecule filename.')
    parser.add_argument('output',
                        help='Output SDF molecule filename.')
    parser.add_argument('-n', '--n-conformers', type=int, default=1,
                        help='Maximum number of conformers to generate for '
                             'each molecule.')
    parser.add_argument('-t', '--rmsd-threshold', type=float, default=0.5,
                        help='RMSD threshold for distinguishin conformers.')
    parser.add_argument('-np', '--n-jobs', type=int, default=1,
                        help='Number of parallel jobs.')
    args = parser.parse_args()
    main()
