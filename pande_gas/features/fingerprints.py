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
    size : int, optional (default 2048)
        Length of generated bit vector.
    chiral : bool, optional (default False)
        Whether to consider chirality in fingerprint generation.
    bonds : bool, optional (default True)
        Whether to consider bond order in fingerprint generation.
    features : bool, optional (default False)
        Whether to use feature information instead of atom information; see
        RDKit docs for more info.
    sparser : bool, optional (default False)
        Whether to return a dict for each molecule containing the sparse
        fingerprint.
    """
    name = 'circular'

    def __init__(self, radius=2, size=2048, chiral=False, bonds=True,
                 features=False, sparse=False):
        self.radius = radius
        self.size = size
        self.chiral = chiral
        self.bonds = bonds
        self.features = features
        self.sparse = sparse

    def _featurize(self, mol):
        """
        Calculate circular fingerprint.

        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        if self.sparse:
            fp = rdMolDescriptors.GetMorganFingerprint(
                mol, self.radius, useChirality=self.chiral,
                useBondTypes=self.bonds, useFeatures=self.features)
            fp = fp.GetNonzeroElements()  # convert to a dict
        else:
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
                mol, self.radius, nBits=self.size, useChirality=self.chiral,
                useBondTypes=self.bonds, useFeatures=self.features)
        return fp
