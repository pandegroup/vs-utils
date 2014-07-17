"""
Test featurize.py.
"""
import cPickle
from cStringIO import StringIO
import gzip
import numpy as np
import os
from rdkit_utils import serial, conformers
import tempfile

from pande_gas.scripts.featurize import main, parse_args


def setup():
    """
    Set up for tests. Writes molecules and targets to files.
    """
    # write mols
    mols = serial.read_mols(StringIO(test_smiles), mol_format='smi')
    mols = [conformers.generate_conformers(mol) for mol in mols]
    _, input_filename = tempfile.mkstemp(suffix='.sdf')
    serial.write_mols_to_file(mols, input_filename)

    # write targets
    targets = [0, 1]
    _, targets_filename = tempfile.mkstemp(suffix='.pkl')
    with open(targets_filename, 'wb') as f:
        cPickle.dump(targets, f, cPickle.HIGHEST_PROTOCOL)

    return input_filename, targets_filename


def cleanup(filenames):
    """
    Delete temporary files.

    Parameters
    ----------
    filenames : list
        Filenames to delete.
    """
    for filename in filenames:
        os.remove(filename)


def test_circular():
    """Test circular fingerprints."""
    input_filename, targets_filename = setup()

    # run script
    _, output_filename = tempfile.mkstemp(suffix='.pkl.gz')
    input_args = [input_filename, '-t', targets_filename, output_filename,
                  'circular', '--size', '2048']
    args = parse_args(input_args)
    main(args.klass, args.input, args.output, args.targets,
         vars(args.featurizer_kwargs))

    # check output file
    with gzip.open(output_filename) as f:
        data = cPickle.load(f)
    assert data['features'].shape == (2, 2048)
    assert data['y'] == [0, 1]
    assert np.array_equal(data['names'], ['aspirin', 'ibuprofen'])

    # cleanup
    cleanup([input_filename, targets_filename, output_filename])


def test_coulomb_matrix():
    """Test Coulomb matrices."""
    input_filename, targets_filename = setup()

    # run script
    _, output_filename = tempfile.mkstemp(suffix='.pkl.gz')
    input_args = [input_filename, '-t', targets_filename, output_filename,
                  'coulomb_matrix', '--max_atoms', '50']
    args = parse_args(input_args)
    main(args.klass, args.input, args.output, args.targets,
         vars(args.featurizer_kwargs))

    # check output file
    with gzip.open(output_filename) as f:
        data = cPickle.load(f)
    assert data['features'].shape == (2, 1, 1275)
    assert data['y'] == [0, 1]
    assert np.array_equal(data['names'], ['aspirin', 'ibuprofen'])

    # cleanup
    cleanup([input_filename, targets_filename, output_filename])


def test_image_features():
    """Test image features."""
    input_filename, targets_filename = setup()

    # run script
    _, output_filename = tempfile.mkstemp(suffix='.pkl.gz')
    input_args = [input_filename, '-t', targets_filename, output_filename,
                  'image', '--size', '16']
    args = parse_args(input_args)
    main(args.klass, args.input, args.output, args.targets,
         vars(args.featurizer_kwargs))

    # check output file
    with gzip.open(output_filename) as f:
        data = cPickle.load(f)
    assert data['features'].shape == (2, 16, 16, 3)
    assert data['y'] == [0, 1]
    assert np.array_equal(data['names'], ['aspirin', 'ibuprofen'])

    # cleanup
    cleanup([input_filename, targets_filename, output_filename])


def test_esp():
    """Test ESP."""
    input_filename, targets_filename = setup()

    # run script
    _, output_filename = tempfile.mkstemp(suffix='.pkl.gz')
    input_args = [input_filename, '-t', targets_filename, output_filename,
                  'esp', '--size', '20']
    args = parse_args(input_args)
    main(args.klass, args.input, args.output, args.targets,
         vars(args.featurizer_kwargs))

    # check output file
    with gzip.open(output_filename) as f:
        data = cPickle.load(f)
    assert data['features'].shape == (2, 1, 61, 61, 61)
    assert data['y'] == [0, 1]
    assert np.array_equal(data['names'], ['aspirin', 'ibuprofen'])

    # cleanup
    cleanup([input_filename, targets_filename, output_filename])

test_smiles = """CC(=O)OC1=CC=CC=C1C(=O)O aspirin
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O ibuprofen
"""
