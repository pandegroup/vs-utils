#!/usr/bin/env python
"""
Generate coulomb matrices for molecules.

See Montavon et al., _New Journal of Physics_ __15__ (2013) 095003.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import cPickle

from rdkit_utils import serial

from pande_gas.features.coulomb_matrices import CoulombMatrix
from pande_gas.utils import h5


def main():
    """Generate coulomb matrices for molecules."""
    mols = serial.read_mols_from_file(args.input)
    featurizer = CoulombMatrix(randomize=args.randomize,
                               n_samples=args.n_samples, seed=args.seed)
    x = featurizer(mols)
    data = {'coulomb_matrices': x}
    if args.labels is not None:
        with open(args.labels) as f:
            y = cPickle.load(f)
        assert len(y) == len(x), 'Labels do not match data.'
        data['y'] = y
    h5.dump(data, args.output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        help='Input SDF molecule filename.')
    parser.add_argument('output',
                        help='Output HDF5 filename.')
    parser.add_argument('-l', '--labels',
                        help='Target labels.')
    parser.add_argument('--no-randomize', dest='randomize',
                        action='store_false',
                        help='Do not randomize Coulomb matrices.')
    parser.add_argument('--n-samples', type=int, default=1,
                        help='Number of random samples.')
    parser.add_argument('--seed', type=int,
                        help='Random seed.')
    args = parser.parse_args()
    main()
