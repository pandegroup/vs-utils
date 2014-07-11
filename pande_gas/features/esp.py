"""
Compute electrostatic potential (ESP) features for molecules.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

from cStringIO import StringIO
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdGeometry, rdMolTransforms

from pande_gas.features import Featurizer
from pande_gas.utils import amber_utils
from pande_gas.utils.pdb_utils import PdbReader


class ESP(Featurizer):
    """
    Calculate electrostatic potential (ESP) features for molecules.

    Parameters
    ----------
    size : float, optional (default 30.)
        Length of each side of the grid, in Angstroms. Used to calculate
        PBSA parameters xmin, xmax, etc.
    resolution : float, optional (default 0.5)
        Space between grid points, in Angstroms. Corresponds to PBSA space
        parameter.
    nb_cutoff : float, optional (default 5.)
        Cutoff distance for van der Waals interactions. Corresponds to PBSA
        cutnb parameter.
    ionic_strength : float, optional (default 150.)
        Ionic strength of the solvent, in mM. Corresponds to PBSA istrng
        parameter.
    """
    conformers = True
    name = 'esp'

    def __init__(self, size=30., resolution=0.5, nb_cutoff=5.,
                 ionic_strength=150.):
        self.size = float(size)
        self.resolution = float(resolution)
        self.nb_cutoff = float(nb_cutoff)
        self.ionic_strength = float(ionic_strength)

    def _featurize(self, mol):
        """
        Calculate electrostatic potential grid.

        Parameters
        ----------
        mol : RDMol
            Molecule.
        """
        return self.calculate_esp(mol)

    def calculate_esp(self, mol):
        """
        Calculate electrostatic potential grid.

        Procedure
        ---------
        1. Prepare molecule by adding hydrogens and canonicalizing the
            orientation and alignment.
        2. Calculate charges and radii with Antechamber.
        3. Calculate electrostatic potential grids with PBSA.

        PBSA requires PQR input, which is similar to PDB with charge and
        radius information added. Neither Antechamber nor PBSA take piped
        input, so we have to use a lot of temporary files and directories.

        Parameters
        ----------
        mol : RDMol
            Molecule.
        """
        mol = self.prepare_molecule(mol)

        # the molecule must have at least one conformer
        assert mol.GetNumConformers()

        # calculate charges and radii
        antechamber = amber_utils.Antechamber()
        charges, radii = antechamber.get_charges_and_radii(mol)

        # get ESP grid for each conformer
        grids = []
        pbsa = amber_utils.PBSA(self.size, self.resolution, self.nb_cutoff,
                                self.ionic_strength)
        for conf in mol.GetConformers():

            # generate a PQR file for this conformer
            pdb = Chem.MolToPDBBlock(mol, confId=conf.GetId())
            reader = PdbReader()
            pqr = reader.pdb_to_pqr(StringIO(pdb), charges, radii)

            # calculate ESP grid
            grid = pbsa.get_esp_grid(pqr)
            grids.append(grid)

        grids = np.asarray(grids)
        return grids

    def prepare_molecule(self, mol):
        """
        Prepare molecule for electrostatic potential calculation.

        * Add hydrogens (with coordinates), since RDKit readers strip
            hydrogens by default.
        * Center and align the molecule.

        Parameters
        ----------
        mol : RDMol
            Molecule.
        """
        mol = Chem.AddHs(mol, addCoords=True)
        center = rdGeometry.Point3D(0, 0, 0)
        for conf in mol.GetConformers():
            rdMolTransforms.CanonicalizeConformer(conf, center=center)
        return mol
