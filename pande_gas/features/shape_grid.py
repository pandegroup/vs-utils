"""
Grid-based shape features.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import numpy as np

from . import Featurizer, MolPreparator
from .gridmol.molecule import GridAtom, GridMol


class ShapeGrid(Featurizer):
    """
    Grid-based shape features.

    Parameters
    ----------
    size : int, optional (default 81)
        Length of each side of the grid, in points.
    resolution : float, optional (default 0.5)
        Space between grid points, in Angstroms.
    hydrogens : bool, optional (default False)
        Whether to consider hydrogens when calculating shape features.
    align : bool, optional (default False)
        Whether to canonicalize the orientation of molecules. This requires
        removal and readdition of hydrogens. This is usually not required
        when working with conformers retrieved from PubChem.
    probe_radius : float, optional (default 1.4)
        Probe radius for determining solvent-accessible surface.
    featurization : str, optional (default 'distance')
        Shape featurization to use. Choose from:
        * 'distance' : distances to molecular surface.
        * 'occupancy' : boolean grid indicating inside/outside of molecule
    distance_to_surface : bool, optional (default True)
        Whether to calculate the distance from each grid point to the molecular
        surface. If False, a boolean grid is returned with set bits
        corresponding to points inside the molecule.
    """
    conformers = True
    name = 'shape'

    def __init__(self, size=81, resolution=0.5, hydrogens=False, align=False,
                 probe_radius=1.4, featurization='distance'):
        self.size = size
        self.resolution = resolution
        self.hydrogens = hydrogens
        self.preparator = MolPreparator(align=align, add_hydrogens=hydrogens)
        self.probe_radius = probe_radius
        self.featurization = featurization

    def _featurize(self, mol):
        """
        Generate shape features for all conformers of a molecule.

        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        mol = self.preparator(mol)
        features = []
        for conf in mol.GetConformers():
            grid_mol = self.embed_mol_in_grid(mol, conf.GetId())
            if self.featurization == 'distance':
                this_features = grid_mol.get_distance()
            elif self.featurization == 'occupancy':
                this_features = grid_mol.get_occupancy()
            else:
                raise NotImplementedError(
                    "Unrecognized featurization '{}'.".format(
                        self.featurization))
            features.append(this_features)
        return features

    def embed_mol_in_grid(self, mol, conf_id):
        """
        Add atoms from a molecule to a GridMol.

        Parameters
        ----------
        mol : RDKit Mol
            RDKit molecule.
        conf_id : int
            RDKit molecule conformer ID.
        """
        shape = tuple(self.size * np.ones(3, dtype=int))
        grid_mol = GridMol(shape, spacing=self.resolution,
                           probe_radius=self.probe_radius)
        conf = mol.GetConformer(conf_id)
        for atom in mol.GetAtoms():
            if not self.hydrogens and atom.GetAtomicNum() == 1:
                continue
            center = list(conf.GetAtomPosition(atom.GetIdx()))
            radius = GridAtom.get_radius_from_atomic_num(atom.GetAtomicNum())
            grid_mol.add_atom(center, radius)
        return grid_mol
