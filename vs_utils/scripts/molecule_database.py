"""
Add molecules to a master database.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse

from rdkit import Chem

from vs_utils.utils.dataset_utils import MoleculeDatabase
from vs_utils.utils.rdkit_utils import serial


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
    parser.add_argument('-d', '--database',
                        help='Existing database to update.')
    parser.add_argument('--stereo-from-3d', action='store_true',
                        help='Assign stereochemistry from 3D coordinates.')
    return parser.parse_args(input_args)


def main(input_filenames, output_filename, database_filename=None,
         assign_stereo_from_3d=False):
    """
    Update or create a molecule database.

    Parameters
    ----------
    input_filenames : list
        Input molecule filename(s).
    output_filename : str
        Output filename.
    database_filename : str, optional
        Existing database to update.
    assign_stereo_from_3d : bool, optional (default False)
        Whether to assign stereochemistry from 3D coordinates.
    """
    database = MoleculeDatabase(assign_stereo_from_3d=assign_stereo_from_3d)
    if database_filename is not None:
        database.load(database_filename)
    initial_size = len(database)
    for filename in input_filenames:
        print filename
        with serial.MolReader().open(filename) as reader:
            for mol in reader:
                try:
                    database.add_mol(mol)
                except ValueError:
                    if mol.HasProp('_Name'):
                        print 'Skipping {}'.format(mol.GetProp('_Name'))
                    else:
                        print 'Skipping {}'.format(
                            Chem.MolToSmiles(mol, isomericSmiles=True))
    final_size = len(database)
    print '{} molecules added to the database'.format(
        final_size - initial_size)
    database.save(output_filename)


if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.output, args.database, args.stereo_from_3d)
