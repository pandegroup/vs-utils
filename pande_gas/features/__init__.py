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
    from .esp import ESP
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
        Whether features are calculated for conformers. If True, the first
        two axes of the feature matrix will index molecules and conformers,
        respectively. If False, only molecule-level features are calculated
        and the feature matrix will not have a separate conformer dimension.
        This is a class attribute because some featurizers take 3D
        conformation into account and others do not, and this is not
        typically an instance-level decision.
    name : str or list
        Name (or names) of this featurizer (used for scripting).
    topo_view : bool (default False)
        Whether the calculated features represent a topological view of the
        data.
    """
    conformers = False
    name = None
    topo_view = False

    def featurize(self, mols, parallel=False, client_kwargs=None,
                  view_flags=None):
        """
        Calculate features for molecules.

        Parameters
        ----------
        mols : iterable
            RDKit Mol objects.
        parallel : bool, optional
            Whether to train subtrainers in parallel using
            IPython.parallel (default False).
        client_kwargs : dict, optional
            Keyword arguments for IPython.parallel Client.
        view_flags : dict, optional
            Flags for IPython.parallel LoadBalancedView.
        """
        if self.conformers and isinstance(mols, types.GeneratorType):
            mols = list(mols)

        if parallel:
            from IPython.parallel import Client

            if client_kwargs is None:
                client_kwargs = {}
            if view_flags is None:
                view_flags = {}
            client = Client(**client_kwargs)
            client.direct_view().use_dill()  # use dill
            view = client.load_balanced_view()
            view.set_flags(**view_flags)
            call = view.map(self._featurize, mols, block=False)
            features = call.get()
        else:
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

    def __call__(self, mols, parallel=False, client_kwargs=None,
                 view_flags=None):
        """
        Calculate features for molecules.

        Parameters
        ----------
        mols : iterable
            RDKit Mol objects.
        parallel : bool, optional
            Whether to train subtrainers in parallel using
            IPython.parallel (default False).
        client_kwargs : dict, optional
            Keyword arguments for IPython.parallel Client.
        view_flags : dict, optional
            Flags for IPython.parallel LoadBalancedView.
        """
        return self.featurize(mols, parallel, client_kwargs, view_flags)

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
        shape = (len(mols), max_confs) + features[0][0].shape
        x = np.ma.masked_all(shape)

        # fill in the container
        for i, (mol, mol_features) in enumerate(zip(mols, features)):
            n_confs = max(mol.GetNumConformers(), 1)
            x[i, :n_confs] = mol_features

        return x
