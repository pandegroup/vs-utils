"""
Feature-based datasets.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import numpy as np

from pylearn2.datasets import DenseDesignMatrix
from pylearn2.space import CompositeSpace
from pylearn2.utils.iteration import FiniteDatasetIterator, safe_izip

from pande_gas.utils import rdkit_utils as rd


class FeaturizerDataset(DenseDesignMatrix):
    """
    Subclass of DenseDesignMatrix that calculates molecular features and
    constructs a dataset.

    Parameters
    ----------
    mols : iterable or str
        Iterable containing RDKit Mol objects or filename containing
        molecules.
    featurizers : list
        FeaturizerModule objects to be applied to molecules.
    y : str, optional
        Key into HDF5 file for dataset targets.
    one_hot : bool, optional (default False)
        If true, convert targets to one-hot encoding.
    mol_iterator : bool, optional (default True)
        Whether to use FeaturizerDatasetIterator as the dataset iterator,
        which indexes the dataset by molecules and not conformers. This is
        required for accurate cross-validation, but must be disabled for
        training without cross-validation.
    kwargs : dict, optional
        Keyword arguments passed to `DenseDesignMatrix`.
    """
    def __init__(self, mols, featurizers, y=None, one_hot=False,
                 mol_iterator=True, **kwargs):
        if isinstance(mols, str):
            mols = rd.read(mols, return_generator=False)

        # all featurizers must have the same setting for topo_view
        topo_view = [featurizer.topo_view for featurizer in featurizers]
        if all(topo_view):
            use_topo_view = True
        else:
            assert all(~np.asarray(topo_view))
            use_topo_view = False

        # if any features are calculated for conformers, we may need to
        # adjust the shape of some featurizers
        conformers = any([featurizer.conformers for featurizer in featurizers])
        if conformers:
            n_confs = [mol.GetNumConformers() for mol in mols]
        else:
            n_confs = None
            mol_iterator = False

        # calculate features
        features = []
        for featurizer in featurizers:
            this_features = featurizer(mols)
            if conformers and not featurizer.conformers:
                this_features = add_conformer_dimension(this_features, n_confs)
            features.append(this_features)

        # construct dataset
        features = np.concatenate(features, axis=-1)
        if conformers and not mol_iterator:
            features = reshape_features(features)
        if use_topo_view:
            x = None
            topo_view = features
        else:
            x = features
            topo_view = None
        super(FeaturizerDataset, self).__init__(X=x, topo_view=topo_view, y=y,
                                                **kwargs)
        if one_hot and y is not None:
            self.convert_to_one_hot()

        # update data_specs to match formatted batch
        if mol_iterator:
            dim = np.prod(features.shape[2:])
            if isinstance(self.data_specs[0], CompositeSpace):
                self.data_specs[0].components[0].dim = dim
            else:
                self.data_specs[0].dim = dim
        self.mol_iterator = mol_iterator

    def iterator(self, *args, **kwargs):
        """
        Get an iterator for this dataset.

        Change the iterator class to FeaturizerDatasetIterator.

        Parameters
        ----------
        WRITEME
        """
        iterator = super(FeaturizerDataset, self).iterator(*args, **kwargs)
        if self.mol_iterator:
            iterator.__class__ = FeaturizerDatasetIterator
        return iterator


class FeaturizerDatasetIterator(FiniteDatasetIterator):
    """
    Dataset iterator that stores data internally indexed by molecules and
    conformers on first two axes, respectively. Selections with __getitem__
    return 2D design matrices, while all other iteration functions index
    the data as it is stored internally.

    Parameters
    ----------
    WRITEME
    """
    def next(self):
        """
        Get the next subset of the dataset during dataset iteration.
        """
        next_index = self._subset_iterator.next()
        rval = []
        for i, (data, fn) in enumerate(
                safe_izip(self._raw_data, self._convert)):
            data = data[next_index]
            if i == 0:
                rval.append(FeaturizerDatasetSelection(data, fn))
            else:
                if fn is not None:
                    data = fn(data)
                rval.append(data)
        rval = tuple(rval)
        if not self._return_tuple and len(rval) == 1:
            rval, = rval
        return rval


class FeaturizerDatasetSelection(object):
    """
    Allow selections on molecules that are then expanded into design
    matrices including conformers. This is useful in cross-validation when
    we need to be able to divide datasets without splitting conformers from
    the same molecule into the training and test data.

    Parameters
    ----------
    data : ndarray
        Molecular features with first two axes indexing molecules and
        conformers, respectively.
    convert : callable, optional
        Callable to be applied to each selection.
    """
    def __init__(self, data, convert=None):
        self.data = data
        self.convert = convert

    def __getitem__(self, item):
        """
        Select molecules and return a design matrix including each set of
        conformer features as a separate example.

        Parameters
        ----------
        item : slice or ndarray
            Molecule selection. Either a slice or a boolean mask.
        """
        rval = self.data[item]
        rval = reshape_features(rval)
        if self.convert is not None:
            rval = self.convert(rval)
        return rval


def add_conformer_dimension(x, n_confs):
    """
    Add conformer dimension to features calculated on molecules.

    Parameters
    ----------
    x : ndarray
        Molecular features.
    n_confs : list
       List containing the number of conformers for each molecule.
    """
    assert all(n_confs)
    features = np.ma.masked_all((x.shape[0], 1) + x.shape[1:])
    for i, n in enumerate(n_confs):
        features[i, :n] = x[i]
    return features


def reshape_features(features):
    """
    Reshape calculated features indexed with molecules and conformers
    as the first two axes, respectively, to match standard 2D or 4D
    representations.

    Parameters
    ----------
    features : ndarray
        Molecular features with first two axes indexing molecules and
        conformers, respectively.
    """
    s = features.shape
    new_shape = (np.prod(s[:2]),) + s[2:]
    rval = np.ma.reshape(features, new_shape)
    return rval
