"""
Basic molecular features.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

from rdkit.Chem import rdMolDescriptors

from pande_gas.features import Featurizer


class MolecularWeight(Featurizer):
    """
    Molecular weight.
    """
    name = ['mw', 'molecular_weight']

    def _featurize(self, mol):
        """
        Calculate molecular weight.

        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        wt = rdMolDescriptors.CalcExactMolWt(mol)
        wt = [wt]
        return wt
