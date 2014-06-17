"""
Topological fingerprints.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

from rdkit.Chem import rdMolDescriptors

from pande_gas.features import Featurizer


class CircularFingerprint(Featurizer):
    """
    Circular (Morgan) fingerprints.

    Parameters
    ----------
    radius : int, optional (default 2)
        Fingerprint radius.
    size : int, optional (default 1024)
        Length of generated bit vector.
    chiral : bool, optional (default False)
        Whether to consider chirality in fingerprint generation.
    bonds : bool, optional (default True)
        Whether to consider bond order in fingerprint generation.
    features : bool, optional (default False)
        Whether to use feature information instead of atom information; see
        RDKit docs for more info.
    """
    name = 'circular'

    def __init__(self, radius=2, size=1024, chiral=False, bonds=False,
                 features=False):
        self.radius = radius
        self.size = size
        self.chiral = chiral
        self.bonds = bonds
        self.features = features

    def _featurize(self, mol):
        """
        Calculate circular fingerprint.

        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol, self.radius, nBits=self.size, useChirality=self.chiral,
            useBondTypes=self.bonds, useFeatures=self.features)
        return fp
