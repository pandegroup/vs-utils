"""
Compute electrostatic potential (ESP) features for molecules.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import numpy as np
import subprocess
import warnings

from rdkit import Chem

from . import Featurizer, MolPreparator
from ..utils import amber_utils
from ..utils.ob_utils import IonizerError


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
    ionize : bool, optional (default True)
        Whether to ionize molecules prior to calculation of ESP.
    pH : float, optional (default 7.4)
        Ionization pH.
    align : bool, optional (default False)
        Whether to canonicalize the orientation of molecules. This requires
        removal and readdition of hydrogens. This is usually not required
        when working with conformers retrieved from PubChem.
    """
    conformers = True
    name = 'esp'

    def __init__(self, size=30., resolution=0.5, nb_cutoff=5.,
                 ionic_strength=150., ionize=True, pH=7.4, align=False):
        self.size = float(size)
        self.resolution = float(resolution)
        self.nb_cutoff = float(nb_cutoff)
        self.ionic_strength = float(ionic_strength)
        self.preparator = MolPreparator(ionize, pH, align, add_hydrogens=True)

    def _featurize(self, mol):
        """
        Calculate electrostatic potential grid.

        Parameters
        ----------
        mol : RDMol
            Molecule.
        """

        # catch ioniziation failures (turn off ionization and try again)
        try:
            prepared_mol = self.preparator(mol)
        except IonizerError:
            if mol.HasProp('_Name'):
                name = mol.GetProp('_Name')
            else:
                name = Chem.MolToSmiles(mol, isomericSmiles=True)
            warnings.warn("Ionization failed for molecule '{}'.".format(name))
            self.preparator.set_ionize(False)
            prepared_mol = self.preparator(mol)

        # catch subprocess failures (turn off ionization and try again)
        try:
            rval = self.calculate_esp(prepared_mol)
        except subprocess.CalledProcessError as e:
            if self.preparator.ionize:
                if mol.HasProp('_Name'):
                    name = mol.GetProp('_Name')
                else:
                    name = Chem.MolToSmiles(mol, isomericSmiles=True)
                warnings.warn("Disabling ionization for molecule '{}'.".format(
                    name))
                self.preparator.set_ionize(False)
                prepared_mol = self.preparator(mol)
                rval = self.calculate_esp(prepared_mol)
            else:
                raise e
        return rval

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

        # the molecule must have at least one conformer
        if mol.GetNumConformers() < 1:
            name = ''
            if mol.HasProp('_Name'):
                name = mol.GetProp('_Name')
            raise AssertionError(
                "Molecule '{}' has zero conformers.".format(name))

        # calculate charges and radii
        antechamber = amber_utils.Antechamber()
        charges, radii = antechamber.get_charges_and_radii(mol)

        # get ESP grid for each conformer
        grids = []
        pbsa = amber_utils.PBSA(self.size, self.resolution, self.nb_cutoff,
                                self.ionic_strength)
        for conf in mol.GetConformers():
            grid, center = pbsa.get_esp_grid(mol, charges, radii,
                                             conf_id=conf.GetId())
            assert center == (0, 0, 0)  # should be centered on the origin
            grids.append(grid)

        grids = np.asarray(grids)
        return grids
