"""
Molecule images.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

from rdkit.Chem import Draw

from pande_gas.features import Featurizer
from pande_gas.utils import image_utils


class MolImage(Featurizer):
    """
    Molecule images.

    Parameters
    ----------
    size : int, optional (default 32)
        Size (in any dimension) of generated images.
    flatten : bool, optional (default False)
        Whether to flatten the pixel array. If False, the features for each
        molecule will be a 3D array.
    """
    name = 'image'

    def __init__(self, size=32, flatten=False):
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
        dim = (self.size, self.size)
        image = Draw.MolToImage(mol, dim, fitImage=True)
        pixels = image_utils.get_pixels(image)
        if self.flatten:
            pixels = pixels.ravel()
        return pixels
