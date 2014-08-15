"""
Molecules on grids.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "3-clause BSD"

import numpy as np
from scipy.spatial.distance import cdist

from . import Grid


class GridMol(Grid):
    """
    Molecule on a grid.

    Keeps track of the distance of each grid point to the molecular surface.
    Points inside the molecule are assigned negative values.

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
    probe_radius : float, optional (default 1.4)
        Probe radius for determining solvent-accessible surface.
    """
    def __init__(self, shape, center=None, spacing=1., dtype=float,
                 probe_radius=1.4):
        super(GridMol, self).__init__(shape, center, spacing, dtype)
        self.probe_radius = probe_radius
        self.atoms = []

    def get_num_atoms(self):
        """
        Get the number of atoms in the molecule.
        """
        return len(self.atoms)

    def add_atom(self, center, radius):
        """
        Add an atom to the molecule.

        Parameters
        ----------
        center : tuple
            Coordinates of atom center.
        radius : float
            Atomic radius.
        """
        atom = GridAtom(self, center, radius)
        self.atoms.append(atom)

    def get_boolean_grid(self):
        """
        Get a boolean grid with set bits corresponding to points inside the
        molecule.
        """
        grid = np.sum([atom.get_grid_mask() for atom in self.atoms])
        grid = np.asarray(grid, dtype=bool)
        return grid

    def get_distance_to_surface(self):
        """
        Get distances to the molecular surface.

        The distance to the molecular surface is the difference between the
        distance to the atom center and the radius of the atom (plus the probe
        radius). The chosen value for each grid point is the minimum absolute
        distance calculated with respect to any atom in the molecule.

        This definition assigns negative values to points within the molecule.
        """
        assert self.get_num_atoms()
        coords = self.get_all_coords().reshape((self.size, self.ndim))
        centers = np.asarray([atom.center for atom in self.atoms], dtype=float)
        radii = np.asarray([atom.radius for atom in self.atoms], dtype=float)
        distances = cdist(coords, centers)
        distances -= (radii + self.probe_radius)

        # take the minimum absolute distance while preserving sign
        min_abs = np.argmin(np.fabs(distances), axis=1)
        distances = distances[np.arange(distances.shape[0]), min_abs]
        distances = distances.reshape(self.shape)
        return distances


class GridAtom(object):
    """
    Atom on a grid.

    Parameters
    ----------
    parent : Molecule
        Parent molecule.
    center : tuple
        Coordinates of atom center.
    radius : float
        Atomic radius.
    """
    def __init__(self, parent, center, radius):
        if not self.atom_is_in_grid(parent, center, radius,
                                    parent.probe_radius):
            raise ValueError('Atom does not fit in the grid.')
        self.parent = parent
        self.center = center
        self.radius = radius

    @staticmethod
    def get_radius_from_atomic_num(atomic_num):
        """
        Look up Bondi van der Waals radius (in Angstroms) for atoms.

        Radii are taken from DOI: 10.1021/j100785a001.

        Parameters
        ----------
        atomic_num : int
            Atomic number.
        """
        radii = {1: 1.2, 6: 1.7, 7: 1.55, 8: 1.52, 9: 1.47, 15: 1.8, 16: 1.8,
                 17: 1.75, 35: 1.85, 53: 1.98}
        rval = radii[atomic_num]
        return rval

    @staticmethod
    def atom_is_in_grid(grid, center, radius, probe_radius=0.):
        """
        Check that this atom is within the parent grid.

        Parameters
        ----------
        parent : Molecule
            Parent molecule.
        center : tuple
            Coordinates of atom center.
        radius : float
            Atomic radius.
        probe_radius : float, optional (default 0.)
            Probe radius for determining solvent-accessible surface.
        """
        try:
            assert grid.coords_in_grid(center)
            for i in xrange(grid.ndim):
                pos_point = np.array(center)
                pos_point[i] += (radius + probe_radius)
                neg_point = np.array(center)
                neg_point[i] -= (radius + probe_radius)
                assert grid.coords_in_grid(pos_point)
                assert grid.coords_in_grid(neg_point)
            return True
        except AssertionError:
            return False

    def get_grid_mask(self):
        """
        Get a boolean mask for grid points inside this atom.
        """
        coords = self.parent.get_all_coords().reshape(
            (self.parent.size, self.parent.ndim))
        distance = cdist(coords, np.atleast_2d(self.center))
        distance = distance.reshape(self.parent.shape)
        mask = distance <= (self.radius + self.parent.probe_radius)
        return mask
