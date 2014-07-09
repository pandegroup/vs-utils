"""
Molecule images.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import os
import subprocess

from rdkit import Chem

from pande_gas.features import Featurizer
from pande_gas.utils import image_utils


class MolImage(Featurizer):
    """
    Molecule images.

    Parameters
    ----------
    size : int, optional (default 32)
        Size (in any dimension) of generated images.
    flatten : bool, optional (default True)
        Whether to flatten the pixel array. If False, the features for each
        molecule will be a 3D array.
    """
    name = 'image'

    def __init__(self, size=32, flatten=True):
        self.size = size
        if not flatten:
            self.topo_view = True
        self.flatten = flatten

    def _featurize(self, mol):
        """
        Generate a 2D depiction of a molecule.

        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        image = self.image_from_smiles(smiles)
        pixels = image_utils.get_pixels(image)
        if self.flatten:
            pixels = pixels.ravel()
        return pixels

    def image_from_smiles(self, smiles):
        """
        Generate a PNG image from a SMILES string using OpenBabel.

        Note that we call OpenBabel using Popen to avoid entanglements with
        the GNU GPL.

        See http://stackoverflow.com/questions/13332268 for details on the
        echo function.

        Parameters
        ----------
        smiles : str
            SMILES string.
        """
        devnull = open(os.devnull, 'w')
        echo = lambda string: subprocess.Popen(['echo', string],
                                               stdout=subprocess.PIPE).stdout
        png_args = ['obabel', '-ican', '-opng', '-xd', '-xC',
                    '-xp {}'.format(self.size)]
        png = subprocess.check_output(png_args, stdin=echo(smiles),
                                      stderr=devnull)
        im = image_utils.load(png)
        return im
