"""
Tests for amber_utils.
"""
from cStringIO import StringIO
import numpy as np
import unittest

from rdkit import Chem
from rdkit_utils import conformers

from vs_utils.utils import amber_utils, pdb_utils


class TestAmberUtils(unittest.TestCase):
    """
    Tests for amber_utils.
    """
    def setUp(self):
        """
        Set up tests.
        """
        smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
        mol = Chem.MolFromSmiles(smiles)
        engine = conformers.ConformerGenerator(max_conformers=1)
        self.mol = engine.generate_conformers(mol)
        assert self.mol.GetNumConformers() > 0


class TestAntechamber(TestAmberUtils):
    """
    Tests for Antechamber.
    """
    def test_antechamber_charges_and_radii(self):
        """
        Test Antechamber charges and radii.
        """
        antechamber = amber_utils.Antechamber()
        charges, radii = antechamber.get_charges_and_radii(self.mol)
        assert charges.size == 21  # 21 atoms: C_9H_8O_4
        assert np.allclose(charges.sum(), 0)  # neutral molecule
        assert radii.size == 21  # 12 atoms
        assert np.count_nonzero(radii > 0)  # no zero radii


class TestPBSA(TestAmberUtils):
    """
    Test PBSA.
    """
    def setUp(self):
        """
        Set up tests.
        """
        super(TestPBSA, self).setUp()

        # get charges and radii
        antechamber = amber_utils.Antechamber()
        self.charges, self.radii = antechamber.get_charges_and_radii(self.mol)

    def test_mol_to_pqr(self):
        """
        Test PBSA.mol_to_pqr.
        """
        pdb = """HEADER    First atom from 4NIP (modified)
ATOM      1  N   GLY A   1       4.168   1.038   6.389  1.00 21.86           N
END
"""
        pqr = """COMPND    First atom from 4NIP (modified)
ATOM 1 N GLY A 1 4.168 1.038 6.389 -0.35 2.3
END
"""
        mol = Chem.MolFromPDBBlock(pdb)
        assert amber_utils.PBSA.mol_to_pqr(mol, [-0.35], [2.3]) == pqr

    def test_pbsa_esp_grid(self):
        """
        Test PBSA.get_esp_grid.
        """
        # calculate ESP grid
        pbsa = amber_utils.PBSA()
        grid, _ = pbsa.get_esp_grid(self.mol, self.charges, self.radii)

        # the grid should be cubic
        size = grid.shape[0]
        assert grid.shape == (size, size, size), grid.shape

        # and not be all zeros
        assert np.count_nonzero(grid)

    def test_pbsa_esp_grid_from_pqr(self):
        """
        Test PBSA.get_esp_grid_from_pqr.
        """
        # write PQR
        reader = pdb_utils.PdbReader()
        pdb = Chem.MolToPDBBlock(self.mol)
        pqr = reader.pdb_to_pqr(StringIO(pdb), self.charges, self.radii)

        # calculate ESP grid
        pbsa = amber_utils.PBSA()
        grid, _ = pbsa.get_esp_grid_from_pqr(pqr)

        # the grid should be cubic
        size = grid.shape[0]
        assert grid.shape == (size, size, size), grid.shape

        # and not be all zeros
        assert np.count_nonzero(grid)
