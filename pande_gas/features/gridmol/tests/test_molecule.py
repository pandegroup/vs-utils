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

    def test_get_occupancy(self):
        """
        Test GridMol.get_occupancy.
        """
        self.mol.add_atom((1, 2, 1), 1.6)
        self.mol.add_atom((1, 1, 1), 1.6)
        occupancy = self.mol.get_occupancy()
        assert occupancy.shape == self.mol.shape
        assert np.count_nonzero(occupancy == 0)
        assert np.count_nonzero(occupancy == 1)
        assert not np.count_nonzero(occupancy > 1)

        # check that most of the grid is empty
        assert np.count_nonzero(occupancy) < 0.2 * self.mol.size

    def test_get_distance(self):
        """
        Test GridMol.get_distance.
        """
        self.mol.add_atom((1, 2, 1), 1.6)
        distances = self.mol.get_distance()

        # confirm that negative values are inside atoms
        mask = self.mol.atoms[0].get_grid_mask()
        assert np.all(distances[mask] <= 0)
        assert np.all(distances[~mask] > 0)

        # check for sane positive distances
        assert np.amax(distances) < max(self.mol.get_real_shape())

        # check that negative distances are not too large
        # since there is only one atom, min should be no larger than the atom
        # radius plus the probe radius
        assert np.fabs(np.amin(distances)) <= (
            self.mol.atoms[0].radius + self.mol.probe_radius)

        # check that most distances are significantly less than max
        threshold = max(self.mol.get_real_shape()) / 2.
        assert np.count_nonzero(np.fabs(distances) < threshold) > (
            0.9 * distances.size)


class TestGridAtom(unittest.TestCase):
    """
    Tests for GridAtom.
    """
    def setUp(self):
        """
        Set up tests.
        """
        self.mol = GridMol((11, 11, 11))
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
