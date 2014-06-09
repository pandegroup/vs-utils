"""
Feature calculations.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import numpy as np


class Featurizer(object):
    """
    Abstract class for calculating a set of features for a molecule.

    Child classes implement the _featurize method for calculating features
    for a single molecule. The feature matrix returned by featurize has a
    shape that is inferred from the shape of the features returned by
    _featurize, and is either indexed by molecules or by molecules and
    conformers depending on the value of the conformers class attribute.

    Class Attributes
    ----------------
    conformers : bool, optional (default False)
        Whether features are calculated for all conformers. If True, the
        first two axes of the feature matrix will index molecules and
        conformers, respectively. If False, features will be calculated
        once for each molecule.
    topo_view : bool (default False)
        Whether the calculated features represent a topological view of the
        data.
    """
    conformers = False
    topo_view = False

    def featurize(self, mols):
        """
        Calculate features for molecules.

        Parameters
        ----------
        mols : iterable
            RDKit Mol objects.
        """
        x = None
        for i, mol in enumerate(mols):
            features = self._featurize(mol)
            if not hasattr(features, 'shape'):
                features = np.asarray(features)
            if x is None:
                x = self.get_container(mols, features.shape)
            if self.conformers:
                x[i, :max(mol.GetNumConformers(), 1)] = features
            else:
                x[i] = features
        return x

    def _featurize(self, mol):
        """
        Calculate features for a single molecule.

        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        raise NotImplementedError('Featurizer is not defined.')

    def __call__(self, mols):
        """
        Calculate features for molecules.

        Parameters
        ----------
        mols : iterable
            RDKit Mol objects.
        """
        return self.featurize(mols)

    def get_container(self, mols, shape):
        """
        Initialize a masked feature container.

        Parameters
        ----------
        mols : iterable
            RDKit Mol objects.
        shape : tuple or int
            Shape or number of features for each conformer.
        """
        if isinstance(shape, int):
            shape = (shape,)
        if self.conformers:
            max_confs = max([mol.GetNumConformers() for mol in mols])
            if not max_confs:
                max_confs = 1
            # max_confs replaces shape[0]
            shape = (len(mols), max_confs) + shape[1:]
        else:
            shape = (len(mols),) + shape
        rval = np.ma.masked_all(shape)
        return rval
