"""
Molecule images.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

from rdkit.Chem import Draw

from vs_utils.features import Featurizer
from vs_utils.utils import image_utils, ob_utils


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
    engine : str, optional (default 'obabel')
        Which engine to use to generate images. Choose from 'obabel' or
        'rdkit'.
    """
    name = 'image'

    def __init__(self, size=32, flatten=False, engine='obabel'):
        self.size = size
        if not flatten:
            self.topo_view = True
        self.flatten = flatten
        self.engine = engine

    def _featurize(self, mol):
        """
        Generate a 2D depiction of a molecule.

        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        dim = (self.size, self.size)
        if self.engine == 'obabel':
            image = ob_utils.MolImage(self.size)(mol)
        elif self.engine == 'rdkit':
            image = Draw.MolToImage(mol, dim, fitImage=True)
            image = image.convert('RGB')  # drop alpha channel
        else:
            raise NotImplementedError(self.engine)
        pixels = image_utils.get_pixels(image)
        if self.flatten:
            pixels = pixels.ravel()
        return pixels
