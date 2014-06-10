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
from pande_gas.utils import images


class MolImage(Featurizer):
    """
    Molecule images.

    Parameters
    ----------
    shape : tuple of ints, optional
        Shape of generated images.
    flatten : bool, optional (default False)
        Whether to flatten the pixel array.
    """
    def __init__(self, shape=None, flatten=True):
        self.shape = shape
        if not flatten:
            self.topo_view = True
        self.flatten = flatten

    def _featurize(self, mol):
        """
        Generate an image for the given molecule and return pixels.

        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        smiles = Chem.MolToSmiles(mol)
        im = self.image_from_smiles(smiles)
        im = images.pad(im, self.shape)
        pixels = images.get_pixels(im)
        if self.flatten:
            pixels = pixels.ravel()
        return pixels

    def image_from_smiles(self, smiles):
        """
        Generate a PNG image from a SMILES string. Uses OpenBabel to
        generate an SVG (which gives uniform scaling) and ImageMagick
        convert to generate a PNG.

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
        svg_args = ['obabel', '-ican', '-osvg', '-d', '-xd']
        svg = subprocess.check_output(svg_args, stdin=echo(smiles),
                                      stderr=devnull)
        png_args = ['convert', '-alpha', 'off', 'svg:-', 'png:-']
        png = subprocess.check_output(png_args, stdin=echo(svg),
                                      stderr=devnull)
        im = images.get_image_from_string(png, max_size=min(self.shape))
        return im
