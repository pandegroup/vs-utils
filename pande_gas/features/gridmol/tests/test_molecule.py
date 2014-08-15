"""
Tests for molecule.py.
"""
import numpy as np
import unittest

from ..molecule import GridAtom, GridMol


class TestGridMol(unittest.TestCase):
    """
    Tests for GridMol.
    """
    def setUp(self):
        """
        Set up tests.
        """
        self.mol = GridMol((11, 11, 11))

    def test_add_atom(self):
        """
        Test GridMol.add_atom.
        """
        self.mol.add_atom((1, 2, 1), 1.6)
        assert len(self.mol.atoms) == 1

    def test_get_boolean_grid(self):
        """
        Test GridMol.get_boolean_grid.
        """
        self.mol.add_atom((1, 2, 1), 1.6)
        self.mol.add_atom((1, 3, 1), 1.6)
        grid = self.mol.get_boolean_grid()
        assert np.count_nonzero(grid == 0)
        assert np.count_nonzero(grid == 1)
        assert not np.count_nonzero(grid > 1)

    def test_get_distance_to_surface(self):
        """
        Test GridMol.update_distance_to_surface.
        """
        self.mol.add_atom((1, 2, 1), 1.6)
        distances = self.mol.get_distance_to_surface()

        # confirm that negative values are inside atoms
        mask = self.mol.atoms[0].get_grid_mask()
        assert np.all(distances[mask] <= 0)
        assert np.all(distances[~mask] > 0)


class TestGridAtom(TestGridMol):
    """
    Tests for GridAtom.
    """
    def setUp(self):
        """
        Set up tests.
        """
        super(TestGridAtom, self).setUp()
        self.atom = GridAtom(self.mol, (1, 2, 1), 1.6)

    def test_atom_is_in_grid(self):
        """
        Test GridAtom.atom_is_in_grid.
        """
        assert GridAtom.atom_is_in_grid(self.mol, self.atom.center,
                                        self.atom.radius,
                                        self.mol.probe_radius)
        assert not GridAtom.atom_is_in_grid(self.mol, (1, 2, 3), 1.6,
                                            self.mol.probe_radius)

    def test_get_grid_mask(self):
        """
        Test GridAtom.get_grid_mask.
        """
        mask = self.atom.get_grid_mask()
        assert np.count_nonzero(mask)

        # add up grid volume and compare to atom volume
        grid_volume = np.count_nonzero(mask) * self.mol.spacing ** 3
        effective_radius = self.atom.radius + self.mol.probe_radius
        atom_volume = 4/3. * np.pi * effective_radius ** 3
        assert np.fabs(grid_volume - atom_volume) < 10
