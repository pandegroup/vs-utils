"""
Test shard_dataset.py.
"""
import shutil
import tempfile
import unittest

from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit_utils import serial

from pande_gas.scripts.shard_dataset import get_filenames, main, write


class TestShardDataset(unittest.TestCase):
    """
    Test shard_dataset.py.
    """
    def setUp(self):
        """
        Set up tests.
        """
        self.reader = serial.MolReader()

        # generate molecules
        smiles = ['CC(=O)OC1=CC=CC=C1C(=O)O', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
                  'CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F']
        names = ['aspirin', 'ibuprofen', 'celecoxib']
        self.mols = []
        for s, n in zip(smiles, names):
            mol = Chem.MolFromSmiles(s)
            mol.SetProp('_Name', n)
            AllChem.Compute2DCoords(mol)
            self.mols.append(mol)

        # write molecules to file
        self.temp_dir = tempfile.mkdtemp()
        writer = serial.MolWriter()
        _, self.filename = tempfile.mkstemp(dir=self.temp_dir,
                                            suffix='.sdf.gz')
        with writer.open(self.filename) as w:
            w.write(self.mols)

    def tearDown(self):
        """
        Clean up tests.
        """
        shutil.rmtree(self.temp_dir)

    def compare_mols(self, filename, ref_slice=None):
        """
        Compare sharded molecules with original molecules.

        Parameters
        ----------
        filename : str
            Filename containing sharded molecules.
        ref_slice : slice, optional
            Slice of self.mols to compare with sharded molecules.
        """
        ref_mols = self.mols
        if ref_slice is not None:
            ref_mols = self.mols[ref_slice]
        mols = list(self.reader.open(filename).get_mols())
        assert len(mols) == len(ref_mols)
        for a, b in zip(mols, ref_mols):
            assert Chem.MolToSmiles(a) == Chem.MolToSmiles(b)
            assert a.GetProp('_Name') == b.GetProp('_Name')

    def test_write(self):
        """
        Test write.
        """
        _, filename = tempfile.mkstemp(dir=self.temp_dir, suffix='.sdf.gz')
        write(self.mols, filename)
        self.compare_mols(filename)

    def test_preserve_properties(self):
        """
        Test preservation of molecule properties when pickling.
        """
        _, filename = tempfile.mkstemp(dir=self.temp_dir, suffix='.pkl.gz')
        write(self.mols, filename)
        self.compare_mols(filename)

    def test_get_filenames(self):
        """
        Test get_filenames.
        """
        filenames = get_filenames('foo', 'bar', index=5)
        for i in xrange(10):
            assert filenames.next() == 'foo-{}.bar'.format(i + 5)

    def test_main(self):
        """
        Test main.
        """
        main(self.filename, output_prefix='foo')
        self.compare_mols('foo-0.pkl.gz')

    def test_leftover(self):
        """
        Test total % chunk_size != 0.
        """
        main(self.filename, chunk_size=2, output_prefix='foo')
        self.compare_mols('foo-0.pkl.gz', slice(2))
        self.compare_mols('foo-1.pkl.gz', slice(2, 3))
