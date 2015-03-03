"""
Tests for molecule_database.py.
"""
import shutil
import tempfile
import unittest

from rdkit import Chem

from rdkit_utils import conformers, serial

from vs_utils.scripts.molecule_database import main, parse_args
from vs_utils.utils.dataset_utils import MoleculeDatabase


class TestMoleculeDatabase(unittest.TestCase):
    """
    Tests for molecule_database.py.
    """
    def setUp(self):
        """
        Set up tests.
        """
        smiles = ['CC(=O)OC1=CC=CC=C1C(=O)O', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
                  'CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F']
        names = ['aspirin', 'ibuprofen', 'celecoxib']
        self.cids = [2244, 3672, 2662]
        self.mols = []
        for s, n in zip(smiles, names):
            mol = Chem.MolFromSmiles(s)
            mol.SetProp('_Name', n)
            self.mols.append(mol)
        self.temp_dir = tempfile.mkdtemp()
        _, self.input_filename = tempfile.mkstemp(dir=self.temp_dir,
                                                  suffix='.smi')
        with open(self.input_filename, 'wb') as f:
            for this_smiles, name in zip(smiles, names):
                f.write('{}\t{}\n'.format(this_smiles, name))
        _, self.output_filename = tempfile.mkstemp(dir=self.temp_dir)

    def tearDown(self):
        """
        Clean up tests.
        """
        shutil.rmtree(self.temp_dir)

    def check_output(self, input_args):
        """
        Run main and examine the resulting database.

        Parameters
        ----------
        args : list
            Command-line arguments.
        """
        args = parse_args(input_args)
        main(args.input, args.output, args.database, args.stereo_from_3d)
        database = MoleculeDatabase()
        database.load(args.output)
        assert len(database) == len(self.mols)
        return database

    def test_defaults(self):
        """
        Test default arguments.
        """
        self.check_output(
            ['-i', self.input_filename, '-o', self.output_filename])

    def test_update(self):
        """
        Test updating an existing database.
        """
        _, database_filename = tempfile.mkstemp(dir=self.temp_dir)
        database = MoleculeDatabase()
        database.add_mol(self.mols[0])
        database.save(database_filename)
        self.check_output(
            ['-i', self.input_filename, '-o', self.output_filename, '-d',
             database_filename])

    def test_assign_stereo_from_3d(self):
        """
        Test --stereo-from-3d.
        """
        # generate conformers for ibuprofen
        engine = conformers.ConformerGenerator()
        mol = engine.generate_conformers(self.mols[1])
        assert mol.GetNumConformers() > 0
        self.mols[1] = mol

        # rewrite input file
        _, self.input_filename = tempfile.mkstemp(dir=self.temp_dir,
                                                  suffix='.sdf')
        with serial.MolWriter().open(self.input_filename) as writer:
            writer.write(self.mols)

        # check for absence of chirality using default arguments
        database = self.check_output(
            ['-i', self.input_filename, '-o', self.output_filename])
        chiral = False
        for smiles in database:
            if '@' in smiles:
                chiral = True
        assert not chiral

        # check for presence of chiraliy using --stereo-from-3d
        database = self.check_output(
            ['-i', self.input_filename, '-o', self.output_filename,
             '--stereo-from-3d'])
        chiral = False
        for smiles in database:
            if '@' in smiles:
                chiral = True
        assert chiral
