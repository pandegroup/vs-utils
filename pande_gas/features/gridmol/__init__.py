"""
N-dimensional grids.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "3-clause BSD"

import numpy as np


class Grid(object):
    """
    N-dimensional grid.

    Parameters
    ----------
    shape : tuple
        Number of grid points in each dimension.
    center : tuple, optional (defaults to the origin)
        Grid center.
    spacing : float, optional (defaults to 1.)
        Space between grid points.
    dtype : numpy dtype, optional (defaults to float)
        Grid data type.
    """
    def __init__(self, shape, center=None, spacing=1., dtype=float):
        self.grid = np.zeros(shape, dtype=dtype)
        self.shape = self.grid.shape
        self.size = self.grid.size
        self.ndim = self.grid.ndim
        if center is None:
            self.center = np.zeros_like(self.shape, dtype=float)
        else:
            self.center = np.asarray(center, dtype=float)
        self.spacing = float(spacing)

    def __getitem__(self, item):
        rval = self.grid[item]
        return rval

    def __setitem__(self, key, value):
        self.grid[key] = value

    def get_real_shape(self):
        """
        Get the shape of the grid in real-space.
        """
        return tuple(np.asarray(self.shape) * self.spacing)

    def grid_point_in_grid(self, grid_point):
        """
        Check whether a grid point is in the grid.

        Parameters
        ----------
        grid_point : tuple
            Grid point.
        """

        # get the maximum indices and see if they exist
        test_point = np.amax(np.atleast_2d(grid_point), axis=0)
        try:
            _ = self.grid[test_point]
            return True
        except IndexError:
            return False

    def coords_in_grid(self, coords):
        """
        Check whether a real-space point is in the grid.

        Parameters
        ----------
        coords : tuple
            Real-space coordinates.
        """
        try:
            self.get_grid_point(coords)
        except IndexError:
            return False
        else:
            return True

    def get_coords(self, grid_points):
        """
        Get real-space coordinates for one or more grid points.

        Procedure
        ---------
        1. Get the real-space coordinates of the grid point relative to the
            grid origin (grid[0, ..., 0]).
        2. Get the real-space translation of the grid center relative to the
            grid origin.
        3. Subtract (2) from (1) to get the real-space coordinates of the grid
            point relative to the grid center rather than the grid origin.
        4. Translate the coordinates from (3) so they are relative to the
            real-space grid center.

        Parameters
        ----------
        grid_points : array_like
            Grid point(s).
        """
        assert np.atleast_2d(grid_points).ndim == 2
        if not self.grid_point_in_grid(grid_points):
            raise IndexError("Invalid grid point '{}'".format(grid_points))
        coords = np.array(grid_points, dtype=float) * self.spacing
        grid_center = (np.asarray(self.shape) - 1) / 2. * self.spacing
        coords -= grid_center  # move from grid_origin to grid_center
        coords += self.center  # move to real-space center
        return coords

    def get_all_coords(self):
        """
        Get real-space coordinates for all grid points.
        """

        # construct an array whose entries correspond to each index of the grid
        # e.g. grid_points[1, 2, 3] -> [1, 2, 3]
        grid_points = np.zeros((self.grid.shape + (self.grid.ndim,)),
                               dtype=int)
        for i, indices in enumerate(np.indices(grid_points.shape[:-1])):
            grid_points[..., i] = indices

        # reshape so grid_points is 2D and each row is a grid point
        grid_points = grid_points.reshape((self.size, self.ndim))

        # get real-space coordinates for all grid points
        # reshape to match shape of the grid
        # e.g. coords[1, 2, 3] -> self.get_coords([1, 2, 3])
        coords = self.get_coords(grid_points).reshape(
            (self.grid.shape + (self.grid.ndim,)))
        return coords

    def get_grid_point(self, coords):
        """
        Get the grid point(s) closest to a set of real-space coordinates.

        Parameters
        ----------
        coords : array_like
            Real-space coordinates.
        """
        coords = np.array(coords, dtype=float)
        coords -= self.center
        center = (np.asarray(self.shape) - 1) / 2. * self.spacing
        coords += center
        coords /= self.spacing
        grid_points = np.asarray(np.rint(coords), dtype=int)
        if not self.grid_point_in_grid(grid_points):
            raise IndexError("Invalid grid point '{}'".format(grid_points))
        return grid_points
