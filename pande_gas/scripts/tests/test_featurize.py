"""
Test featurize.py.
"""
import cPickle
import gzip
import joblib
import numpy as np
from rdkit_utils import conformers, serial
import shutil
import tempfile
import unittest

from rdkit import Chem
from rdkit.Chem import AllChem

from pande_gas.scripts.featurize import main, parse_args


class TestFeaturize(unittest.TestCase):
    """
    Test featurize.py.
    """
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

    def test_pickle(self):
        """
        Save features to a pickle.
        """
        # run script
        _, output_filename = tempfile.mkstemp(suffix='.pkl',
                                              dir=self.temp_dir)
        input_args = [self.input_filename, '-t', self.targets_filename,
                      output_filename, 'circular', '--size', '2048']
        args = parse_args(input_args)
        main(args.klass, args.input, args.output, args.targets,
             vars(args.featurizer_kwargs))

        # check output file
        with open(output_filename) as f:
            data = cPickle.load(f)
        assert data['features'].shape == (2, 2048)
        assert data['y'] == [0, 1]
        assert np.array_equal(data['names'], ['aspirin', 'ibuprofen'])

    def test_compressed_pickle(self):
        """
        Save features to a compressed pickle.
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

    def test_joblib(self):
        """
        Save features using joblib.dump.
        """
        # run script
        _, output_filename = tempfile.mkstemp(suffix='.joblib',
                                              dir=self.temp_dir)
        input_args = [self.input_filename, '-t', self.targets_filename,
                      output_filename, 'circular', '--size', '2048']
        args = parse_args(input_args)
        main(args.klass, args.input, args.output, args.targets,
             vars(args.featurizer_kwargs))

        # check output file
        data = joblib.load(output_filename)
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

    def test_shape_grid(self):
        """
        Test ShapeGrid.
        """

        # run script
        _, output_filename = tempfile.mkstemp(suffix='.pkl.gz')
        input_args = [self.input_filename, '-t', self.targets_filename,
                      output_filename, 'shape', '--size', '40']
        args = parse_args(input_args)
        main(args.klass, args.input, args.output, args.targets,
             vars(args.featurizer_kwargs))

        # check output file
        with gzip.open(output_filename) as f:
            data = cPickle.load(f)
        assert data['features'].shape == (2, 1, 40, 40, 40)
        assert data['y'] == [0, 1]
        assert np.array_equal(data['names'], ['aspirin', 'ibuprofen'])

    def test_scaffolds(self):
        """
        Test scaffold generation.
        """
        # run script
        _, output_filename = tempfile.mkstemp(suffix='.pkl',
                                              dir=self.temp_dir)
        input_args = [self.input_filename, '-t', self.targets_filename,
                      output_filename, 'circular']
        args = parse_args(input_args)
        main(args.klass, args.input, args.output, args.targets,
             vars(args.featurizer_kwargs))

        # check output file
        with open(output_filename) as f:
            data = cPickle.load(f)

        assert Chem.MolFromSmiles(data['scaffolds'][0]).GetNumAtoms() == 6
        assert Chem.MolFromSmiles(data['scaffolds'][1]).GetNumAtoms() == 6

    def test_chiral_scaffolds(self):
        """
        Test chiral scaffold generation.
        """

        # romosetron
        mol = Chem.MolFromSmiles(
            'CN1C=C(C2=CC=CC=C21)C(=O)[C@@H]3CCC4=C(C3)NC=N4')
        AllChem.Compute2DCoords(mol)
        self.mols[1] = mol

        # write mols
        _, self.input_filename = tempfile.mkstemp(suffix='.sdf',
                                                  dir=self.temp_dir)
        writer = serial.MolWriter()
        writer.open(self.input_filename)
        writer.write(self.mols)
        writer.close()

        # run script w/o chiral scaffolds
        _, output_filename = tempfile.mkstemp(suffix='.pkl',
                                              dir=self.temp_dir)
        input_args = [self.input_filename, '-t', self.targets_filename,
                      output_filename, 'circular']
        args = parse_args(input_args)
        main(args.klass, args.input, args.output, args.targets,
             vars(args.featurizer_kwargs),
             chiral_scaffolds=args.chiral_scaffolds)

        # get achiral scaffold
        with open(output_filename) as f:
            data = cPickle.load(f)

        achiral_scaffold = data['scaffolds'][1]

        # run script w/ chiral scaffolds
        input_args = [self.input_filename, '-t', self.targets_filename,
                      output_filename, '--chiral-scaffolds', 'circular']
        args = parse_args(input_args)
        main(args.klass, args.input, args.output, args.targets,
             vars(args.featurizer_kwargs),
             chiral_scaffolds=args.chiral_scaffolds)

        # get chiral scaffold
        with open(output_filename) as f:
            data = cPickle.load(f)

        chiral_scaffold = data['scaffolds'][1]

        assert achiral_scaffold != chiral_scaffold

    def test_collate_mols1(self):
        """
        Test collate_mols where molecules are pruned.
        """

        # write targets
        targets = {'names': ['ibuprofen'], 'y': [0]}
        with open(self.targets_filename, 'wb') as f:
            cPickle.dump(targets, f, cPickle.HIGHEST_PROTOCOL)

        # run script
        _, output_filename = tempfile.mkstemp(suffix='.pkl',
                                              dir=self.temp_dir)
        input_args = [self.input_filename, '-t', self.targets_filename,
                      output_filename, 'circular']
        args = parse_args(input_args)
        main(args.klass, args.input, args.output, args.targets,
             vars(args.featurizer_kwargs))

        # check output file
        with open(output_filename) as f:
            data = cPickle.load(f)

        assert np.array_equal(data['names'], targets['names'])
        assert np.array_equal(data['y'], targets['y'])
        assert data['features'].shape[0] == 1

    def test_collate_mols2(self):
        """
        Test collate_mols where targets are pruned.
        """

        # write targets
        targets = {'names': ['aspirin', 'ibuprofen'], 'y': [0, 1]}
        with open(self.targets_filename, 'wb') as f:
            cPickle.dump(targets, f, cPickle.HIGHEST_PROTOCOL)

        # write mols
        writer = serial.MolWriter()
        writer.open(self.input_filename)
        writer.write([self.mols[0]])
        writer.close()

        # run script
        _, output_filename = tempfile.mkstemp(suffix='.pkl',
                                              dir=self.temp_dir)
        input_args = [self.input_filename, '-t', self.targets_filename,
                      output_filename, 'circular']
        args = parse_args(input_args)
        main(args.klass, args.input, args.output, args.targets,
             vars(args.featurizer_kwargs))

        # check output file
        with open(output_filename) as f:
            data = cPickle.load(f)

        assert np.array_equal(data['names'], ['aspirin'])
        assert np.array_equal(data['y'], [0])
        assert data['features'].shape[0] == 1

    def test_collate_mols3(self):
        """
        Test collate_mols where targets are in a different order than
        molecules.
        """

        # write targets
        targets = {'names': ['ibuprofen', 'aspirin'], 'y': [1, 0]}
        with open(self.targets_filename, 'wb') as f:
            cPickle.dump(targets, f, cPickle.HIGHEST_PROTOCOL)

        # run script
        _, output_filename = tempfile.mkstemp(suffix='.pkl',
                                              dir=self.temp_dir)
        input_args = [self.input_filename, '-t', self.targets_filename,
                      output_filename, 'circular']
        args = parse_args(input_args)
        main(args.klass, args.input, args.output, args.targets,
             vars(args.featurizer_kwargs))

        # check output file
        with open(output_filename) as f:
            data = cPickle.load(f)

        sort = np.argsort(targets['names'])  # names will be sorted
        assert np.array_equal(data['names'],
                              np.asarray(targets['names'])[sort])
        assert np.array_equal(data['y'], np.asarray(targets['y'])[sort])
        assert data['features'].shape[0] == 2
