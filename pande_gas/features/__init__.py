"""
Feature calculations.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import numpy as np
import types


def get_featurizers():
    """Compile a dict mapping strings to featurizer classes."""
    from .basic import MolecularWeight
    from .coulomb_matrices import CoulombMatrix
    from .fingerprints import CircularFingerprint
    from .images import MolImage

    featurizers = {}
    for klass in Featurizer.__subclasses__():
        assert klass.name is not None, (klass.__name__ +
                                        " 'name' attribute is None.")
        if isinstance(klass.name, list):
            for name in klass.name:
                assert name not in featurizers
                featurizers[name] = klass
        else:
            assert klass.name not in featurizers
            featurizers[klass.name] = klass
    return featurizers


def resolve_featurizer(name):
    """
    Resolve featurizer class from a string.

    Parameters
    ----------
    name : str
        Featurizer name.
    """
    return get_featurizers()[name]


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
    name : str or list
        Name (or names) of this featurizer (used for scripting).
    topo_view : bool (default False)
        Whether the calculated features represent a topological view of the
        data.
    """
    conformers = False
    name = None
    topo_view = False

    def featurize(self, mols):
        """
        Calculate features for molecules.

        Parameters
        ----------
        mols : iterable
            RDKit Mol objects.
        """
        if self.conformers and isinstance(mols, types.GeneratorType):
            mols = list(mols)
        features = [self._featurize(mol) for mol in mols]
        if self.conformers:
            features = self.conformer_container(mols, features)
        else:
            features = np.asarray(features)
        return features

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

    def conformer_container(self, mols, features):
        """
        Put features into a container with an extra dimension for
        conformers. Conformer indices that are not used are masked.

        For example, if mols contains 3 molecules with 1, 2, 5 conformers,
        respectively, then the final container will have 3 entries on its
        first axis and 5 entries on its second axis. The remaining axes
        correspond to feature dimensions.

        Parameters
        ----------
        mols : iterable
            RDKit Mol objects.
        features : list
            Features calculated for molecule conformers. Each element
            corresponds to the features for a molecule and should be an
            ndarray with conformers on the first axis.
        """

        # get the maximum number of conformers
        max_confs = max([mol.GetNumConformers() for mol in mols])
        if not max_confs:
            max_confs = 1

        # construct the new container
        # - first axis = # mols
        # - second axis = max # conformers
        # - remaining axes = determined by feature shape
        shape = (len(mols), max_confs) + features[0].shape[1:]
        x = np.ma.masked_all(shape)

        # fill in the container
        for i, (mol, mol_features) in enumerate(zip(mols, features)):
            n_confs = max(mol.GetNumConformers(), 1)
            x[i, :n_confs] = mol_features
            
        return x
