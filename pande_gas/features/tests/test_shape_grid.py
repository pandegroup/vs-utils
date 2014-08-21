"""
Tests for grid-based shape features.
"""
import unittest

from rdkit import Chem
from rdkit_utils import conformers

from ..shape_grid import ShapeGrid


class TestShapeGrid(unittest.TestCase):
    """
    Tests for ShapeGrid featurizer.
    """
    def setUp(self):
        """
        Set up tests.
        """
        smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
        mol = Chem.MolFromSmiles(smiles)
        engine = conformers.ConformerGenerator(max_conformers=3)
        self.mols = [engine.generate_conformers(mol)]
        for mol in self.mols:
            assert mol.GetNumConformers() > 0
        self.max_confs = max([mol.GetNumConformers() for mol in self.mols])
        self.engine = ShapeGrid()

    def test_distance(self):
        """
        Test ShapeGrid with distances featurization.
        """
        self.engine = ShapeGrid(featurization='distance')
        features = self.engine(self.mols)
        assert features.shape == (len(self.mols), self.max_confs, 81, 81, 81)

    def test_occupancy(self):
        """
        Test ShapeGrid with occupancy featurization.
        """
        self.engine = ShapeGrid(featurization='occupancy')
        features = self.engine(self.mols)
        assert features.shape == (len(self.mols), self.max_confs, 81, 81, 81)

    def test_embed_mol_in_grid(self):
        """
        Test ShapeGrid.embed_mol_in_grid.
        """
        grid_mol = self.engine.embed_mol_in_grid(self.mols[0], conf_id=0)
        num_atoms = Chem.RemoveHs(self.mols[0]).GetNumAtoms()
        assert num_atoms < self.mols[0].GetNumAtoms()
        assert grid_mol.get_num_atoms() == num_atoms

    def test_embed_mol_in_grid_hydrogens(self):
        """
        Test ShapeGrid.embed_mol_in_grid with hydrogens=True.
        """
        self.engine = ShapeGrid(hydrogens=True)
        grid_mol = self.engine.embed_mol_in_grid(self.mols[0], conf_id=0)
        num_atoms = self.mols[0].GetNumAtoms()
        assert grid_mol.get_num_atoms() == num_atoms
