"""
Tests for grid.py.
"""
import numpy as np
import unittest

from .. import Grid


class TestGrid(unittest.TestCase):
    """
    Tests for Grid.
    """
    def setUp(self):
        """
        Set up tests.
        """
        self.grid = Grid((11, 11, 11))

    def test_grid_point_in_grid(self):
        """
        Test Grid.grid_point_in_grid.
        """
        assert self.grid.grid_point_in_grid((1, 2, 3))
        assert not self.grid.grid_point_in_grid((11, 2, 3))

    def test_coords_in_grid(self):
        """
        Test Grid.coords_in_grid.
        """
        assert self.grid.coords_in_grid((-1, -1, -1))
        assert not self.grid.coords_in_grid((-1, 6, -1))

    def test_get_coords(self):
        """
        Test Grid.get_coords.
        """
        coords = self.grid.get_coords((1, 2, 3))
        assert np.array_equal(coords, [-4, -3, -2])

    def test_get_all_coords(self):
        """
        Test Grid.get_all_coords.
        """
        coords = self.grid.get_all_coords()
        assert np.array_equal(coords.shape,
                              (self.grid.shape + (self.grid.ndim,)))

    def test_get_grid_point(self):
        """
        Test Grid.get_grid_point.
        """
        grid_point = self.grid.get_grid_point((-1.2, -1.3, -1.4))
        assert np.array_equal(grid_point, [4, 4, 4])

    def test_spacing(self):
        """
        Test an arbitrary spacing.
        """
        self.grid = Grid((11, 11, 11), spacing=0.3)
        assert np.array_equal(self.grid.get_real_shape(), [3.3, 3.3, 3.3])
        grid_point = self.grid.get_grid_point((-1.2, -1.3, -1.4))
        assert np.array_equal(grid_point, [1, 1, 0]), grid_point
