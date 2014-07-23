"""
Open Babel utilities.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import subprocess

from rdkit import Chem

from pande_gas.utils import image_utils


class Ionizer(object):
    """
    Calculate atomic formal charges at the given pH.

    Parameters
    ----------
    pH : float, optional (default 7.4)
        pH at which to calculate formal charges.
    """
    def __init__(self, pH=7.4):
        self.pH = pH

    def __call__(self, mol):
        """
        Ionize a molecule.

        Parameters
        ----------
        mol : RDMol
            Molecule.
        """
        return self.ionize(mol)

    def ionize(self, mol):
        """
        Ionize a molecule.

        Parameters
        ----------
        mol : RDMol
            Molecule.
        """
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        args = ['obabel', '-i', 'can', '-o', 'can', '-p', str(self.pH)]
        p = subprocess.Popen(args, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        ionized_smiles, _ = p.communicate(smiles)
        mol = Chem.MolFromSmiles(ionized_smiles)
        return mol


class MolImage(object):
    """
    Generate 2D depictions of molecules.

    Parameters
    ----------
    size : int, optional (default 32)
        Size (in any direction) of generated images.
    """
    def __init__(self, size=32):
        self.size = size

    def __call__(self, mol):
        """
        Generate a PNG image from a SMILES string.

        Parameters
        ----------
        mol : RDMol
            Molecule.
        size : int, optional (default 32)
            Size (in any direction) of generated image.
        """
        return self.depict(mol)

    def depict(self, mol):
        """
        Generate a PNG image from a SMILES string.

        Parameters
        ----------
        mol : RDMol
            Molecule.
        size : int, optional (default 32)
            Size (in any direction) of generated image.
        """
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        args = ['obabel', '-i', 'can', '-o', 'png', '-xd', '-xC',
                '-xp {}'.format(self.size)]
        p = subprocess.Popen(args, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        png, _ = p.communicate(smiles)
        im = image_utils.load(png)
        return im
