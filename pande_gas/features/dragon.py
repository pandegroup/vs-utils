"""
Dragon descriptors.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

from pande_gas.features import Featurizer
from pande_gas.utils.dragon_utils import Dragon


class DragonDescriptors(Featurizer):
    """
    Calculate Dragon descriptors.

    Parameters
    ----------
    assign_stereo_from_3d : bool, optional (default False)
        Assign stereochemistry from 3D coordinates. This will overwrite any
        existing stereochemistry information on molecules.
    """
    name = 'dragon'

    def __init__(self, assign_stereo_from_3d=False):
        self.engine = Dragon(assign_stereo_from_3d=assign_stereo_from_3d)

    def _featurize(self, mol):
        """
        Calculate Dragon descriptors.

        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        return self.engine.get_descriptors(mol)
