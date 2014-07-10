"""
Test featurize.py.
"""
import cPickle
from cStringIO import StringIO
import gzip
import numpy as np
import os
from rdkit_utils import serial
import tempfile

from pande_gas.scripts.featurize import main, parse_args


def test_circular():
    """Test circular fingerprints."""

    # write mols
    mols = serial.read_mols(StringIO(test_smiles), mol_format='smi')
    _, input_filename = tempfile.mkstemp(suffix='.sdf')
    serial.write_mols_to_file(mols, input_filename)

    # write targets
    targets = [0, 1]
    _, targets_filename = tempfile.mkstemp(suffix='.pkl')
    with open(targets_filename, 'wb') as f:
        cPickle.dump(targets, f, cPickle.HIGHEST_PROTOCOL)

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
    assert data['y'] == targets
    assert np.array_equal(data['names'], ['aspirin', 'ibuprofen'])

    # cleanup
    for filename in [input_filename, targets_filename, output_filename]:
        os.remove(filename)

test_smiles = """CC(=O)OC1=CC=CC=C1C(=O)O aspirin
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O ibuprofen
"""
