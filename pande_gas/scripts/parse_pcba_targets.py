#!/usr/bin/env python
"""
Parse assay data from PubChem BioAssay (PCBA).
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import cPickle
import gzip
import numpy as np
import pandas as pd


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
    subparsers = parser.add_subparsers(dest='subcommand')
    subparsers.add_parser('classification')
    regression = subparsers.add_parser('regression')
    regression.add_argument('columns', nargs='+',
                            help='Column indices to include.')
    return parser.parse_args(input_args)


def read_pcba_data(filename):
    """
    Read PCBA data.

    Parameters
    ----------
    filename : str
        PCBA data filename.
    """
    if filename.endswith('.gz'):
        with gzip.open(filename) as f:
            df = pd.read_csv(f)
    else:
        df = pd.read_csv(filename)

    # remove duplicates
    df = df.drop_duplicates('PUBCHEM_CID')

    return df


def read_map(filename):
    """
    Read ID->SMILES map.

    Parameters
    ----------
    filename : str
        ID->SMILES map filename.
    """
    if filename.endswith('.gz'):
        f = gzip.open(filename)
    else:
        f = open(filename)
    id_map = cPickle.load(f)
    f.close()
    return id_map


def map_smiles(df, id_map):
    """
    Get SMILES for PubChem CIDs.

    Parameters
    ----------
    df : DataFrame
        DataFrame containing PCBA data.
    id_map : dict
        ID->SMILES map.
    """
    smiles = []
    indices = []
    for i, cid in enumerate(df.PUBCHEM_CID):
        name = 'CID{}'.format(cid)  # no bare IDs in map
        try:
            smiles.append(id_map[name])
            indices.append(i)
        except KeyError:
            pass
    return np.asarray(smiles), np.asarray(indices)


def save(smiles, targets, filename):
    """
    Save SMILES and targets to disk.

    Parameters
    ----------
    smiles : array_like
        SMILES for molecules.
    targets : array_like
        Targets for molecules.
    filename : str
        Output filename.
    """
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'wb')
    else:
        f = open(filename, 'wb')
    cPickle.dump({'smiles': smiles, 'targets': targets}, f,
                 cPickle.HIGHEST_PROTOCOL)
    f.close()


def classification_main(input_filename, map_filename, output_filename):
    """
    Get classification targets.

    Parameters
    ----------
    input_filename : str
        PCBA data filename.
    map_filename : str
        ID->SMILES map filename.
    output_filename : str
        Output filename.
    """
    df = read_pcba_data(input_filename)
    id_map = read_map(map_filename)
    smiles, indices = map_smiles(df, id_map)
    targets = np.asarray(df.PUBCHEM_ACTIVITY_OUTCOME == 'Active')[indices]
    save(smiles, targets, output_filename)


def regression_main(input_filename, map_filename, output_filename, cols):
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
    cols : list
        Column indices to include.
    """
    cols = np.asarray(cols, dtype=int)
    df = read_pcba_data(input_filename)
    print "Extracting data from the following columns:"
    for i in cols:
        print '\t', df.columns[i]
    id_map = read_map(map_filename)
    smiles, indices = map_smiles(df, id_map)
    targets = np.zeros((df.shape[0], len(cols)), dtype=float)
    for i, idx in enumerate(cols):
        targets[:, i] = df[df.columns[idx]]
    targets = targets[indices]
    save(smiles, targets, output_filename)

if __name__ == '__main__':
    args = parse_args()
    if args.subcommand == 'classification':
        classification_main(args.input, args.map, args.output)
    elif args.subcommand == 'regression':
        regression_main(args.input, args.map, args.output, args.columns)
