"""
Test featurize.py.
"""
import joblib
import numpy as np
import shutil
import tempfile
import unittest

from rdkit import Chem
from rdkit.Chem import AllChem

from vs_utils.scripts.featurize import main, parse_args
from vs_utils.utils import read_csv_features, read_pickle, write_pickle
from vs_utils.utils.rdkit_utils import conformers, serial


class TestFeaturize(unittest.TestCase):
    """
    Test featurize.py.
    """
    def setUp(self):
        """
        Set up for tests. Writes molecules and targets to files.
        """
        self.temp_dir = tempfile.mkdtemp()
        smiles = ['CC(=O)OC1=CC=CC=C1C(=O)O',
                  'C[C@@H](C1=CC=C(C=C1)CC(C)C)C(=O)O']
        self.names = ['aspirin', 'ibuprofen']
        engine = conformers.ConformerGenerator(max_conformers=1)
        self.mols = []
        self.smiles = []  # use RDKit-generated SMILES
        for i in xrange(len(smiles)):
            mol = Chem.MolFromSmiles(smiles[i])
            mol.SetProp('_Name', self.names[i])
            self.mols.append(engine.generate_conformers(mol))
            self.smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True,
                                                canonical=True))

        # write mols
        _, self.input_filename = tempfile.mkstemp(suffix='.sdf',
                                                  dir=self.temp_dir)
        writer = serial.MolWriter()
        writer.open(self.input_filename)
        writer.write(self.mols)
        writer.close()

        # write targets
        self.targets = [0, 1]
        _, self.targets_filename = tempfile.mkstemp(suffix='.pkl',
                                                    dir=self.temp_dir)
        write_pickle(self.targets, self.targets_filename)

    def tearDown(self):
        """
        Delete temporary files.

        Parameters
        ----------
        filenames : list
            Filenames to delete.
        """
        shutil.rmtree(self.temp_dir)

    def check_output(self, featurize_args, shape, targets=None, names=None,
                     smiles=None, output_suffix='.pkl'):
        """
        Check features shape, targets, and names.

        Parameters
        ----------
        featurize_args : list
            Featurizer-specific arguments for script.
        filename : str
            Output filename.
        shape : tuple
            Expected shape of features.
        targets : list, optional
            Expected targets. Defaults to self.targets.
        names : list, optional
            Expected names. Defaults to self.names.
        smiles : list, optional
            Expected SMILES. Defaults to self.smiles.
        output_suffix : str, optional (default '.pkl')
            Suffix for output files.
        """

        # generate command-line arguments
        _, output_filename = tempfile.mkstemp(suffix=output_suffix,
                                              dir=self.temp_dir)
        input_args = [self.input_filename, '-t', self.targets_filename,
                      output_filename, '--smiles'] + featurize_args

        # run script
        args = parse_args(input_args)
        main(args.klass, args.input, args.output, target_filename=args.targets,
             featurizer_kwargs=vars(args.featurizer_kwargs),
             include_smiles=args.include_smiles, scaffolds=args.scaffolds,
             chiral_scaffolds=args.chiral_scaffolds,
             mol_id_prefix=args.mol_prefix)

        # read output file
        if output_filename.endswith('.joblib'):
            data = joblib.load(output_filename)
        elif (output_filename.endswith('.csv')
                or output_filename.endswith('.csv.gz')):
            data = read_csv_features(output_filename)
        else:
            data = read_pickle(output_filename)

        # check values
        if targets is None:
            targets = self.targets
        if names is None:
            names = self.names
        if smiles is None:
            smiles = self.smiles
        assert len(data) == shape[0]
        if len(shape) > 1:
            assert np.asarray(data.ix[0, 'features']).shape == shape[1:]
        assert np.array_equal(data['y'], targets), data['y']
        assert np.array_equal(data['mol_id'], names), data['mol_id']
        assert np.array_equal(data['smiles'], smiles), data['smiles']

        # return output in case anything else needs to be checked
        return data

    def test_pickle(self):
        """
        Save features to a pickle.
        """
        self.check_output(['circular'], (2, 2048))

    def test_compressed_pickle(self):
        """
        Save features to a compressed pickle.
        """
        self.check_output(['circular'], (2, 2048), output_suffix='.pkl.gz')

    def test_joblib(self):
        """
        Save features using joblib.dump.
        """
        self.check_output(['circular'], (2, 2048), output_suffix='.joblib')

    def test_csv(self):
        """
        Save features to csv.
        """
        self.check_output(['circular'], (2, 2048), output_suffix='.csv')

    def test_csv_gz(self):
        """
        Save features to csv.gz.
        """
        self.check_output(['circular'], (2, 2048), output_suffix='.csv.gz')

    def test_circular(self):
        """
        Test circular fingerprints.
        """
        self.check_output(['circular', '--size', '512'], (2, 512))

    def test_sparse_circular(self):
        """
        Test sparse circular fingerprints.
        """
        data = self.check_output(['circular', '--sparse'], (2,))
        for value in data['features']:
            assert isinstance(value, dict)
            assert len(value)

    def test_coulomb_matrix(self):
        """
        Test Coulomb matrices.
        """
        self.check_output(['coulomb_matrix', '--max_atoms', '50'],
                          (2, 1, 1275))

    def test_image_features(self):
        """
        Test image features.
        """
        self.check_output(['image', '--size', '16'], (2, 16, 16, 3))

    def test_esp(self):
        """
        Test ESP.
        """
        self.check_output(['esp', '--size', '20'], (2, 1, 61, 61, 61))

    def test_shape_grid(self):
        """
        Test ShapeGrid.
        """
        self.check_output(['shape', '--size', '40'], (2, 1, 40, 40, 40))

    def test_mw(self):
        """
        Test calculation of molecular weight.
        """
        self.check_output(['mw'], (2, 1))

    def test_descriptors(self):
        """
        Test calculation of RDKit descriptors.
        """
        self.check_output(['descriptors'], (2, 196))

    def test_mol_id_prefix(self):
        """
        Test that names are prepended with mol_id_prefix.
        """
        prefix = 'CID'
        for i, name in enumerate(self.names):
          self.names[i] = prefix + name
        self.check_output(['--mol-prefix', prefix, 'circular', '--size', '512'],
                          (2, 512))

    def test_scaffold(self):
        """
        Test scaffold featurizer.
        """
        self.check_output(['scaffold'], (2,))

    def test_scaffolds(self):
        """
        Test scaffold generation.
        """
        data = self.check_output(['--scaffolds', 'circular'], (2, 2048))
        assert Chem.MolFromSmiles(data['scaffolds'][0]).GetNumAtoms() == 6
        assert Chem.MolFromSmiles(data['scaffolds'][1]).GetNumAtoms() == 6

    def test_chiral_scaffolds(self):
        """
        Test chiral scaffold generation.
        """

        # romosetron
        mol = Chem.MolFromSmiles(
            'CN1C=C(C2=CC=CC=C21)C(=O)[C@@H]3CCC4=C(C3)NC=N4')
        mol.SetProp('_Name', 'romosetron')
        AllChem.Compute2DCoords(mol)
        self.mols[1] = mol
        self.names[1] = 'romosetron'
        self.smiles[1] = Chem.MolToSmiles(mol, isomericSmiles=True)

        # write mols
        _, self.input_filename = tempfile.mkstemp(suffix='.sdf',
                                                  dir=self.temp_dir)
        writer = serial.MolWriter()
        writer.open(self.input_filename)
        writer.write(self.mols)
        writer.close()

        # run script w/o chiral scaffolds
        data = self.check_output(['--scaffolds', 'circular'], (2, 2048))
        achiral_scaffold = data['scaffolds'][1]

        # run script w/ chiral scaffolds
        data = self.check_output(['--scaffolds', '--chiral-scaffolds',
                                  'circular'], (2, 2048))
        chiral_scaffold = data['scaffolds'][1]

        assert achiral_scaffold != chiral_scaffold

    def test_collate_mols1(self):
        """
        Test collate_mols where molecules are pruned.
        """

        # write targets
        targets = {'names': ['ibuprofen'], 'y': [0]}
        write_pickle(targets, self.targets_filename)

        # run script
        self.check_output(['circular'], (1, 2048), targets=targets['y'],
                          names=targets['names'], smiles=[self.smiles[1]])

    def test_collate_mols2(self):
        """
        Test collate_mols where targets are pruned.
        """

        # write targets
        targets = {'names': ['aspirin', 'ibuprofen'], 'y': [0, 1]}
        write_pickle(targets, self.targets_filename)

        # write mols
        writer = serial.MolWriter()
        writer.open(self.input_filename)
        writer.write([self.mols[0]])
        writer.close()

        # run script
        self.check_output(['circular'], (1, 2048), targets=[0],
                          names=['aspirin'], smiles=[self.smiles[0]])

    def test_collate_mols3(self):
        """
        Test collate_mols where targets are in a different order than
        molecules.
        """

        # write targets
        targets = {'names': ['ibuprofen', 'aspirin'], 'y': [1, 0]}
        write_pickle(targets, self.targets_filename)

        # run script
        self.check_output(['circular'], (2, 2048))
