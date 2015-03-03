"""
Generate molecular scaffolds.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

from vs_utils.features import Featurizer
from vs_utils.utils import ScaffoldGenerator


class Scaffold(Featurizer):
    """
    Molecular scaffolds.

    Parameters
    ----------
    Parameters
    ----------
    include_chirality : : bool, optional (default False)
        Include chirality in scaffolds.
    """
    name = 'scaffold'

    def __init__(self, include_chirality=False):
        self.include_chirality = include_chirality
        self.engine = ScaffoldGenerator(include_chirality=include_chirality)

    def _featurize(self, mol):
        """
        Generate a 2D depiction of a molecule.

        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        return self.engine.get_scaffold(mol)
