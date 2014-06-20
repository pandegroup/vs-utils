"""
Molecule images.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import numpy as np
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
    max_size : int, optional (default 32)
        Maximum size (in any dimension) of generated images.
    flatten : bool, optional (default True)
        Whether to flatten the pixel array. If False, the features for each
        molecule will be a 3D array.
    """
    name = 'image'

    def __init__(self, max_size=32, flatten=True):
        self.max_size = max_size
        if not flatten:
            self.topo_view = True
        self.flatten = flatten

    def featurize(self, mols):
        """
        Generate images for molecules and extract pixels. This class can't
        use _featurize because pad_size needs to be calculated from the
        ensemble of images for mols.

        Parameters
        ----------
        mols : iterable
            RDKit Mol objects.
        """
        smiles = [Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                  for mol in mols]
        images = map(self.image_from_smiles, smiles)
        pad_size = max([max(image.size) for image in images])
        padded = [image_utils.pad(image, (pad_size, pad_size))
                  for image in images]
        scaled = [image_utils.downscale(image, self.max_size)
                  for image in padded]
        pixels = map(image_utils.get_pixels, scaled)
        if self.flatten:
            pixels = [p.ravel() for p in pixels]
        pixels = np.asarray(pixels)
        return pixels

    @staticmethod
    def image_from_smiles(smiles):
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
        im = image_utils.load(png)
        return im
