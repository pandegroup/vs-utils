"""
Construct target files for datasets with active/decoy labels.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import numpy as np

from rdkit_utils import serial

from pande_gas.utils import write_pickle
from pande_gas.utils.dataset_utils import MoleculeDatabase


def parse_args(input_args=None):
    """
    Parse command-line arguments.

    Parameters
    ----------
    input_args : list, optional
        Input arguments. If not provided, defaults to sys.argv[1:].
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--actives', required=1,
                        help='Active molecule filename.')
    parser.add_argument('-d', '--decoys', required=1,
                        help='Decoy molecule filename.')
    parser.add_argument('-o', '--output', required=1,
                        help='Output filename.')
    parser.add_argument('--stereo-from-3d', action='store_true',
                        help='Assign stereochemistry from 3D coordinates.')
    return parser.parse_args(input_args)


def main(active_filename, decoy_filename, output_filename,
         assign_stereo_from_3d=False):
    """
    Construct target files for datasets with active/decoy labels.

    Parameters
    ----------
    active_filename : str
        Active molecule filename.
    decoy_filename : str
        Decoy molecule filename.
    output_filename : str
        Output filename.
    assign_stereo_from_3d : bool, optional (default False)
        Assign stereochemistry from 3D coordinates.
    """
    active_smiles = get_smiles(active_filename, assign_stereo_from_3d)
    decoy_smiles = get_smiles(decoy_filename, assign_stereo_from_3d)
    targets = np.concatenate((np.ones(len(active_smiles), dtype=int),
                              np.zeros(len(decoy_smiles), dtype=int)))
    smiles = np.concatenate((active_smiles, decoy_smiles))
    write_pickle({'smiles': smiles, 'targets': targets}, output_filename)


def get_smiles(filename, assign_stereo_from_3d=False):
    """
    Get SMILES for molecules.

    Parameters
    ----------
    filename : str
        Input molecule filename.
    assign_stereo_from_3d : bool, optional (default False)
        Assign stereochemistry from 3D coordinates.
    """
    database = MoleculeDatabase(assign_stereo_from_3d=assign_stereo_from_3d)
    with serial.MolReader().open(filename) as reader:
        for mol in reader:
            database.add_mol(mol)
    return list(database.smiles)

if __name__ == '__main__':
    args = parse_args()
    main(args.actives, args.decoys, args.output, args.stereo_from_3d)
