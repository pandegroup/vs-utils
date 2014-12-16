"""
Tests for classfication_targets.py.
"""
import shutil
import tempfile
import unittest

from rdkit import Chem

from rdkit_utils import conformers, serial

from pande_gas.scripts.datasets.classification_targets import main, parse_args
from pande_gas.utils import read_pickle, SmilesGenerator


class TestClassificationTargets(unittest.TestCase):
    """
    Tests for classification_targets.py.
    """
    def setUp(self):
        """
        Set up tests.
        """
        smiles = ['CC(=O)OC1=CC=CC=C1C(=O)O', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
                  'CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F']
        names = ['aspirin', 'ibuprofen', 'celecoxib']
        self.y = [0, 1, 0]
        self.mols = []
        for s, n in zip(smiles, names):
            mol = Chem.MolFromSmiles(s)
            mol.SetProp('_Name', n)
            self.mols.append(mol)

        # write active and decoy files
        self.temp_dir = tempfile.mkdtemp()
        _, self.active_filename = tempfile.mkstemp(dir=self.temp_dir,
                                                   suffix='.smi')
        _, self.decoy_filename = tempfile.mkstemp(dir=self.temp_dir,
                                                  suffix='.smi')
        active = open(self.active_filename, 'wb')
        decoy = open(self.decoy_filename, 'wb')
        for this_smiles, name, y in zip(smiles, names, self.y):
            data = '{}\t{}\n'.format(this_smiles, name)
            if y:
                active.write(data)
            else:
                decoy.write(data)
        active.close()
        decoy.close()
        _, self.output_filename = tempfile.mkstemp(dir=self.temp_dir)

        # get SMILES
        self.engine = SmilesGenerator()
        self.smiles = [self.engine.get_smiles(mol) for mol in self.mols]

    def tearDown(self):
        """
        Clean up tests.
        """
        shutil.rmtree(self.temp_dir)

    def check_output(self, input_args):
        """
        Check main output.

        Parameters
        ----------
        input_args : list
            Command-line arguments.
        """
        args = parse_args(input_args)
        main(args.actives, args.decoys, args.output, args.stereo_from_3d)
        data = read_pickle(self.output_filename)
        for smiles, target in zip(data['smiles'], data['targets']):
            assert smiles in self.smiles
            assert target == self.y[self.smiles.index(smiles)]
        return data['smiles'], data['targets']

    def test_defaults(self):
        """
        Test main with default parameters.
        """
        args = ['-a', self.active_filename, '-d', self.decoy_filename, '-o',
                self.output_filename]
        self.check_output(args)

    def test_stereo_to_3d(self):
        """
        Test main with --stereo-to-3d.
        """
        # generate conformers for ibuprofen
        engine = conformers.ConformerGenerator()
        self.mols[1] = engine.generate_conformers(self.mols[1])
        assert self.mols[1].GetNumConformers() > 0

        # rewrite actives file with 3D coordinates
        _, self.active_filename = tempfile.mkstemp(dir=self.temp_dir,
                                                   suffix='.sdf')
        with serial.MolWriter().open(self.active_filename) as writer:
            for mol, y in zip(self.mols, self.y):
                if y:
                    writer.write([self.mols[1]])

        # check for absence of chirality using default arguments
        smiles, targets = self.check_output(
            ['-a', self.active_filename, '-d', self.decoy_filename, '-o',
             self.output_filename])
        chiral = False
        for this_smiles in smiles:
            if '@' in this_smiles:
                chiral = True
        assert not chiral

        # update reference SMILES
        self.engine = SmilesGenerator(assign_stereo_from_3d=True)
        self.smiles[1] = self.engine.get_smiles(self.mols[1])

        # check for presence of chiraliy using --stereo-from-3d
        smiles, targets = self.check_output(
            ['-a', self.active_filename, '-d', self.decoy_filename, '-o',
             self.output_filename, '--stereo-from-3d'])
        chiral = False
        for this_smiles in smiles:
            if '@' in this_smiles:
                chiral = True
        assert chiral
