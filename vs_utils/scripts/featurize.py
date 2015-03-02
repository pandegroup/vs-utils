#!/usr/bin/env python
"""
Featurize molecules and save features to disk. Featurizers are exposed as
subcommands, with __init__ arguments as subcommand arguments.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import argparse
import inspect
import joblib
import numpy as np

from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit_utils import serial

from vs_utils.features import get_featurizers
from vs_utils.utils import (read_pickle, ScaffoldGenerator, SmilesGenerator,
                             write_pickle)
from vs_utils.utils.parallel_utils import LocalCluster


def parse_args(input_args=None):
    """
    Parse command-line arguments. Each featurizer class is a subcommand
    whose arguments are stored in args.featurizer_kwargs. The featurizer
    class is stored in args.klass.

    Parameters
    ----------
    input_args : list, optional
        Input arguments. If not provided, defaults to sys.argv[1:].
    """
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument('input',
                        help='Input molecules.')
    parser.add_argument('--chiral-scaffolds', action='store_true',
                        help='Whether to include chirality in scaffolds.')
    parser.add_argument('-t', '--targets',
                        help='Molecule targets.')
    parser.add_argument('--smiles-hydrogens', action='store_true')
    parser.add_argument('--scaffolds', action='store_true',
                        help='Calculate molecule scaffolds.')
    parser.add_argument('--names', action='store_true',
                        help='Include molecule names.')
    parser.add_argument('-p', '--parallel', action='store_true',
                        help='Whether to use IPython.parallel.')
    parser.add_argument('-id', '--cluster-id',
                        help='IPython.parallel cluster ID.')
    parser.add_argument('-np', '--n-engines', type=int,
                        help='Start a local IPython.parallel cluster with ' +
                             'this many engines.')
    parser.add_argument('output',
                        help='Output filename (.joblib, .pkl, or .pkl.gz).')
    parser.add_argument('-c', '--compression-level', type=int, default=3,
                        help='Compression level (0-9) to use with ' +
                             'joblib.dump.')

    # featurizer subcommands
    featurizers = get_featurizers()
    subparsers = parser.add_subparsers(title='featurizers')
    for name, klass in featurizers.items():
        command = subparsers.add_parser(name, help=klass.__doc__,
                                        formatter_class=HelpFormatter,
                                        epilog=klass.__doc__)
        command.set_defaults(klass=klass)
        try:
            args, _, _, defaults = inspect.getargspec(klass.__init__)
        except TypeError:
            args = []
        for i, arg in enumerate(args):
            if i == 0 and arg == 'self':
                continue
            kwargs = {}
            try:
                kwargs['default'] = defaults[i-len(args)]
                if kwargs['default'] is not None:
                    kwargs['type'] = type(kwargs['default'])
            except IndexError:
                kwargs['required'] = True
            if 'type' in kwargs and kwargs['type'] == bool:
                if kwargs['default']:
                    command.add_argument('--no-{}'.format(arg), dest=arg,
                                         action='store_false')
                else:
                    command.add_argument('--{}'.format(arg),
                                         action='store_true')
            else:
                command.add_argument('--{}'.format(arg), **kwargs)
    args = argparse.Namespace()
    args.featurizer_kwargs = parser.parse_args(input_args)
    for arg in ['input', 'output', 'klass', 'targets', 'parallel',
                'cluster_id', 'n_engines', 'compression_level',
                'smiles_hydrogens', 'names', 'scaffolds', 'chiral_scaffolds']:
        setattr(args, arg, getattr(args.featurizer_kwargs, arg))
        delattr(args.featurizer_kwargs, arg)
    return args


class HelpFormatter(argparse.RawTextHelpFormatter):
    """
    Argparse help formatter with better indenting.

    Parameters
    ----------
    WRITEME
    """
    def __init__(self, prog, indent_increment=2, max_help_position=8,
                 width=None):
        super(HelpFormatter, self).__init__(prog, indent_increment,
                                            max_help_position, width)


def main(featurizer_class, input_filename, output_filename,
         target_filename=None, featurizer_kwargs=None, parallel=False,
         client_kwargs=None, view_flags=None, compression_level=3,
         smiles_hydrogens=False, names=False, scaffolds=False,
         chiral_scaffolds=False):
    """
    Featurize molecules in input_filename using the given featurizer.

    Parameters
    ----------
    featurizer_class : Featurizer
        Featurizer class.
    input_filename : str
        Filename containing molecules to be featurized.
    output_filename : str
        Output filename. Should end with .pkl or .pkl.gz.
    target_filename : str, optional
        Pickle containing target values. Should either be array_like or a dict
        containing 'names' and 'y' keys, corresponding to molecule names and
        target values.
    featurizer_kwargs : dict, optional
        Keyword arguments passed to featurizer.
    parallel : bool, optional
        Whether to train subtrainers in parallel using IPython.parallel
        (default False).
    client_kwargs : dict, optional
        Keyword arguments for IPython.parallel Client.
    view_flags : dict, optional
        Flags for IPython.parallel LoadBalancedView.
    compression_level : int, optional (default 3)
        Compression level (0-9) to use with joblib.dump.
    smiles_hydrogens : bool, optional (default False)
        Whether to keep hydrogens when generating SMILES.
    names : bool, optional (default False)
        Whether to include molecule names in output.
    scaffolds : bool, optional (default False)
        Whether to include scaffolds in output.
    chiral_scaffods : bool, optional (default False)
        Whether to include chirality in scaffolds.
    """
    mols, mol_names = read_mols(input_filename)

    # get targets
    data = {}
    if target_filename is not None:
        targets = read_pickle(target_filename)
        if isinstance(targets, dict):
            mol_indices, target_indices = collate_mols(
                mols, mol_names, targets['y'], targets['names'])
            mols = mols[mol_indices]
            mol_names = mol_names[mol_indices]
            targets = np.asarray(targets['y'])[target_indices]
        else:
            assert len(targets) == len(mols)
        data['y'] = targets

    # featurize molecules
    print "Featurizing molecules..."
    if featurizer_kwargs is None:
        featurizer_kwargs = {}
    featurizer = featurizer_class(**featurizer_kwargs)
    features = featurizer.featurize(mols, parallel, client_kwargs, view_flags)

    # fill in data container
    print "Saving results..."
    data['features'] = features

    # calculate SMILES
    smiles = SmilesGenerator(remove_hydrogens=(not smiles_hydrogens))
    data['smiles'] = np.asarray([smiles.get_smiles(mol) for mol in mols])

    # sanity checks
    assert data['features'].shape[0] == len(mols), (
        "Features do not match molecules.")
    assert data['smiles'].shape[0] == len(mols), (
        "SMILES do not match molecules.")

    # names, scaffolds, args
    if names:
        data['names'] = mol_names
    if scaffolds:
        data['scaffolds'] = get_scaffolds(mols, chiral_scaffolds)
    data['args'] = {'featurizer_class': featurizer_class.__name__,
                    'input_filename': input_filename,
                    'target_filename': target_filename,
                    'featurizer_kwargs': featurizer_kwargs,
                    'chiral_scaffolds': chiral_scaffolds}

    # write output file
    write_output_file(data, output_filename, compression_level)


def collate_mols(mols, mol_names, targets, target_names):
    """
    Prune and reorder mols to match targets.

    Parameters
    ----------
    mols : array_like
        Molecules.
    mol_names : array_like
        Molecule names.
    targets : array_like
        Targets.
    target_names : array_like
        Molecule names corresponding to targets.

    Returns
    -------
    which_mols : array_like
        Indices of molecules to keep, ordered to match targets.
    keep_targets : array_like
        Targets corresponding to selected molecules. This could differ from
        the input targets if some of the targets do not have a corresponding
        molecule (this often happens when a 3D structure cannot be generated
        for a molecule that has target data).
    """
    print "Collating molecules and targets..."
    assert len(mols) == len(mol_names) and len(targets) == len(target_names)

    # make sure dtypes match for names so comparisons will work properly
    target_names = np.asarray(target_names)
    mol_names = np.asarray(mol_names).astype(target_names.dtype)

    # sanity checks
    if np.unique(mol_names).size != mol_names.size:
        raise ValueError('Molecule names must be unique.')
    if np.unique(target_names).size != target_names.size:
        raise ValueError('Molecule names (for targets) must be unique.')

    # get intersection of mol_names and target_names
    shared_names = np.intersect1d(mol_names, target_names)

    # get indices to select those molecules from mols and targets
    mol_indices = np.zeros_like(shared_names, dtype=int)
    target_indices = np.zeros_like(shared_names, dtype=int)
    for i, name in enumerate(shared_names):
        mol_indices[i] = np.where(mol_names == name)[0][0]
        target_indices[i] = np.where(target_names == name)[0][0]
    return mol_indices, target_indices


def read_mols(input_filename):
    """
    Read molecules from an input file and extract names.

    Parameters
    ----------
    input_filename : str
        Filename containing molecules.
    """
    print "Reading molecules..."
    reader = serial.MolReader()
    reader.open(input_filename)
    mols = []
    names = []
    for mol in reader.get_mols():
        mols.append(mol)
        if mol.HasProp('_Name'):
            names.append(mol.GetProp('_Name'))
        else:
            names.append(None)
    reader.close()
    mols = np.asarray(mols)
    names = np.asarray(names)
    return mols, names


def get_scaffolds(mols, include_chirality=False):
    """
    Get Murcko scaffolds for molecules.

    Murcko scaffolds are described in DOI: 10.1021/jm9602928. They are
    essentially that part of the molecule consisting of rings and the linker
    atoms between them.

    Parameters
    ----------
    mols : array_like
        Molecules.
    include_chirality : bool, optional (default False)
        Whether to include chirality in scaffolds.
    """
    print "Generating molecule scaffolds..."
    engine = ScaffoldGenerator(include_chirality=include_chirality)
    scaffolds = []
    for mol in mols:
        scaffolds.append(engine.get_scaffold(mol))
    scaffolds = np.asarray(scaffolds)
    return scaffolds


def write_output_file(data, output_filename, compression_level=3):
    """
    Pickle output data, possibly to a compressed file.

    Parameters
    ----------
    data : object
        Object to pickle in output file.
    output_filename : str
        Output filename. Should end with .joblib, .pkl, or .pkl.gz.
    compression_level : int, optional (default 3)
        Compression level (0-9) to use with joblib.dump.
    """
    if output_filename.endswith('.pkl') or output_filename.endswith('.pkl.gz'):
        write_pickle(data, output_filename)
    elif output_filename.endswith('.joblib'):
        joblib.dump(data, output_filename, compress=compression_level)
    else:
        raise NotImplementedError('Unrecognized output file extension.')

if __name__ == '__main__':
    args = parse_args()

    # start a cluster
    if args.n_engines is not None:
        assert args.cluster_id is None, ('Cluster ID should not be should ' +
                                         'not be specified when starting a' +
                                         'new cluster.')
        cluster = LocalCluster(args.n_engines)
        args.parallel = True
        args.cluster_id = cluster.cluster_id

    # cluster flags
    if args.cluster_id is not None:
        client_kwargs = {'cluster_id': args.cluster_id}
    else:
        client_kwargs = None
    view_flags = {'retries': 1}

    # run main function
    main(args.klass, args.input, args.output, args.targets,
         vars(args.featurizer_kwargs), args.parallel, client_kwargs,
         view_flags, args.compression_level, args.smiles_hydrogens, args.names,
         args.scaffolds, args.chiral_scaffolds)
