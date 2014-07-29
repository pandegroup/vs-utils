"""
Test featurize.py.
"""
import cPickle
import gzip
import numpy as np
from rdkit_utils import conformers, serial
import shutil
import tempfile
import unittest

from rdkit import Chem

from pande_gas.scripts.featurize import main, parse_args


class TestFeaturize(unittest.TestCase):
    def setUp(self):
        """
        Set up for tests. Writes molecules and targets to files.
        """
        self.temp_dir = tempfile.mkdtemp()
        smiles = ['CC(=O)OC1=CC=CC=C1C(=O)O', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O']
        names = ['aspirin', 'ibuprofen']
        engine = conformers.ConformerGenerator(max_conformers=1)
        self.mols = []
        for i in xrange(len(smiles)):
            mol = Chem.MolFromSmiles(smiles[i])
            mol.SetProp('_Name', names[i])
            self.mols.append(engine.generate_conformers(mol))

        # write mols
        _, self.input_filename = tempfile.mkstemp(suffix='.sdf',
                                                  dir=self.temp_dir)
        writer = serial.MolWriter()
        writer.open(self.input_filename)
        writer.write(self.mols)
        writer.close()

        # write targets
        targets = [0, 1]
        _, self.targets_filename = tempfile.mkstemp(suffix='.pkl',
                                                    dir=self.temp_dir)
        with open(self.targets_filename, 'wb') as f:
            cPickle.dump(targets, f, cPickle.HIGHEST_PROTOCOL)

    def tearDown(self):
        """
        Delete temporary files.

        Parameters
        ----------
        filenames : list
            Filenames to delete.
        """
        shutil.rmtree(self.temp_dir)

    def test_circular(self):
        """
        Test circular fingerprints.
        """

        # run script
        _, output_filename = tempfile.mkstemp(suffix='.pkl.gz',
                                              dir=self.temp_dir)
        input_args = [self.input_filename, '-t', self.targets_filename,
                      output_filename, 'circular', '--size', '2048']
        args = parse_args(input_args)
        main(args.klass, args.input, args.output, args.targets,
             vars(args.featurizer_kwargs))

        # check output file
        with gzip.open(output_filename) as f:
            data = cPickle.load(f)
        assert data['features'].shape == (2, 2048)
        assert data['y'] == [0, 1]
        assert np.array_equal(data['names'], ['aspirin', 'ibuprofen'])

    def test_coulomb_matrix(self):
        """
        Test Coulomb matrices.
        """

        # run script
        _, output_filename = tempfile.mkstemp(suffix='.pkl.gz')
        input_args = [self.input_filename, '-t', self.targets_filename,
                      output_filename, 'coulomb_matrix', '--max_atoms', '50']
        args = parse_args(input_args)
        main(args.klass, args.input, args.output, args.targets,
             vars(args.featurizer_kwargs))

        # check output file
        with gzip.open(output_filename) as f:
            data = cPickle.load(f)
        assert data['features'].shape == (2, 1, 1275)
        assert data['y'] == [0, 1]
        assert np.array_equal(data['names'], ['aspirin', 'ibuprofen'])

    def test_image_features(self):
        """
        Test image features.
        """

        # run script
        _, output_filename = tempfile.mkstemp(suffix='.pkl.gz')
        input_args = [self.input_filename, '-t', self.targets_filename,
                      output_filename, 'image', '--size', '16']
        args = parse_args(input_args)
        main(args.klass, args.input, args.output, args.targets,
             vars(args.featurizer_kwargs))

        # check output file
        with gzip.open(output_filename) as f:
            data = cPickle.load(f)
        assert data['features'].shape == (2, 16, 16, 3)
        assert data['y'] == [0, 1]
        assert np.array_equal(data['names'], ['aspirin', 'ibuprofen'])

    def test_esp(self):
        """
        Test ESP.
        """

        # run script
        _, output_filename = tempfile.mkstemp(suffix='.pkl.gz')
        input_args = [self.input_filename, '-t', self.targets_filename,
                      output_filename, 'esp', '--size', '20']
        args = parse_args(input_args)
        main(args.klass, args.input, args.output, args.targets,
             vars(args.featurizer_kwargs))

        # check output file
        with gzip.open(output_filename) as f:
            data = cPickle.load(f)
        assert data['features'].shape == (2, 1, 61, 61, 61)
        assert data['y'] == [0, 1]
        assert np.array_equal(data['names'], ['aspirin', 'ibuprofen'])
