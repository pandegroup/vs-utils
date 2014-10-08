"""
Miscellaneous utility functions.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import cPickle
import gzip
import numpy as np
import os

from rdkit import Chem
from rdkit_utils import PicklableMol, serial


def read_pickle(filename):
    """
    Read pickled data from (possibly gzipped) files.

    Parameters
    ----------
    filename : str
        Filename.
    """
    if filename.endswith('.gz'):
        f = gzip.open(filename)
    else:
        f = open(filename)
    data = cPickle.load(f)
    f.close()
    return data


def write_pickle(data, filename, protocol=cPickle.HIGHEST_PROTOCOL):
    """
    Write data to a (possibly gzipped) pickle.

    Parameters
    ----------
    data : object
        Object to pickle.
    filename : str
        Filename.
    protocol : int, optional (default cPickle.HIGHEST_PROTOCOL)
        Pickle protocol.
    """
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'wb')
    else:
        f = open(filename, 'wb')
    cPickle.dump(data, f, protocol)
    f.close()


class DatasetSharder(object):
    """
    Split a dataset into chunks.

    Parameters
    ----------
    filename : str, optional
        Input filename. One of filename or mols must be provided.
    mols : iterable, optional
        Molecules to shard. One of filename or mols must be provided.
    shard_size : int, optional (default 1000)
        Number of molecules per shard.
    write_shards : bool, optional (default True)
        Whether to automatically write shards to disk.
    prefix : str, optional
        Prefix for output files.
    flavor : str, optional (default 'pkl.gz')
        Output molecule format used as the extension for shard filenames.
    start_index : int, optional (default 0)
        Starting index for shard filenames.
    """
    def __init__(self, filename=None, mols=None, shard_size=1000,
                 write_shards=True, prefix=None, flavor='pkl.gz',
                 start_index=0):
        if filename is None and mols is None:
            raise ValueError('One of filename or mols must be provided.')
        self.filename = filename
        self.mols = mols
        self.shard_size = shard_size
        self.write_shards = write_shards
        if self.filename is not None and prefix is None:
            prefix = self._guess_prefix()
        if write_shards and prefix is None:
            raise ValueError('One of filename or prefix must be provided ' +
                             'when writing shards.')
        self.prefix = prefix
        self.flavor = flavor
        self.index = start_index
        self.writer = serial.MolWriter()

    def _guess_prefix(self):
        """
        Get the prefix from a filename.

        Takes everything in the basename before the first period. For example,
        the prefix for '../foo.bar.gz' is 'foo'.
        """
        return os.path.basename(self.filename).split('.')[0]

    def _next_filename(self):
        """
        Generate the next shard filename.
        """
        if self.prefix is None:
            raise ValueError('Prefix must be provided when writing shards.')
        filename = '{}-{}.{}'.format(self.prefix, self.index, self.flavor)
        self.index += 1
        return filename

    def read_mols_from_file(self):
        """
        Read molecules from a file.
        """
        with serial.MolReader().open(self.filename) as reader:
            for mol in reader.get_mols():
                yield mol

    def shard(self):
        """
        Split a dataset into chunks.

        If self.write_shards is False, a shard generator is returned. Each
        shard is an ndarray with dtype=object, which gives convenient access
        to ndarray operations (like fancy indexing) for downstream
        applications.
        """
        if self.write_shards:
            for shard in self._shard():
                self.write_shard(shard)
        else:
            return self._shard()

    def _shard(self):
        """
        Split a dataset into chunks.
        """
        if self.mols is None:
            self.mols = self.read_mols_from_file()
        shard = []
        for mol in self.mols:
            shard.append(mol)
            if len(shard) >= self.shard_size:
                yield np.asarray(shard)  # ndarray with dtype=object
                shard = []
        if len(shard):
            yield np.asarray(shard)

    def __iter__(self):
        """
        Iterate through shards.
        """
        return self._shard()

    def write_shard(self, mols):
        """
        Write molecules to the next shard file.

        Molecules are converted to PicklableMols prior to writing to preserve
        properties such as molecule names.

        Parameters
        ----------
        mols : array_like
            Molecules.
        """
        mols = [PicklableMol(mol) for mol in mols]  # preserve properties
        filename = self._next_filename()
        with self.writer.open(filename) as f:
            f.write(mols)


def pad_array(x, shape, fill=0, both=False):
    """
    Pad an array with a fill value.

    Parameters
    ----------
    x : ndarray
        Matrix.
    shape : tuple or int
        Desired shape. If int, all dimensions are padded to that size.
    fill : object, optional (default 0)
        Fill value.
    both : bool, optional (default False)
        Whether to split the padding on both sides of each axis. If False,
        padding is applied to the end of each axis.
    """
    x = np.asarray(x)
    if not isinstance(shape, tuple):
        shape = tuple(shape for _ in xrange(x.ndim))
    pad = []
    for i in xrange(x.ndim):
        diff = shape[i] - x.shape[i]
        assert diff >= 0
        if both:
            a, b = divmod(diff, 2)
            b += a
            pad.append((a, b))
        else:
            pad.append((0, diff))
    pad = tuple(pad)
    x = np.pad(x, pad, mode='constant', constant_values=fill)
    return x


class SmilesGenerator(object):
    """
    Generate SMILES strings for molecules.

    Parameters
    ----------
    remove_hydrogens : bool, optional (default True)
        Whether to remove hydrogens prior to generating SMILES.
    assign_stereo_from_3d : bool, optional (default False)
        Whether to assign stereochemistry from 3D coordinates. This will
        overwrite any existing stereochemistry information on molecules.
    """
    def __init__(self, remove_hydrogens=True, assign_stereo_from_3d=False):
        self.remove_hydrogens = remove_hydrogens
        self.assign_stereo_from_3d = assign_stereo_from_3d

    def get_smiles(self, mol):
        """
        Map a molecule name to its corresponding SMILES string.

        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        if self.assign_stereo_from_3d:  # do this before removing hydrogens
            Chem.AssignAtomChiralTagsFromStructure(mol)
        if self.remove_hydrogens:
            mol = Chem.RemoveHs(mol)  # creates a copy
        return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)

    def get_unique_smiles(self, mols):
        """
        Get unique SMILES for a set of molecules.

        Parameters
        ----------
        mols : iterable
            Molecules.
        """
        return np.unique([self.get_smiles(mol) for mol in mols])


class SmilesMap(object):
    """
    Map compound names to SMILES.

    Parameters
    ----------
    prefix : str, optional
        Prefix to prepend to IDs.
    allow_duplicates : bool, optional (default True)
        Whether to allow duplicate SMILES.
    kwargs : dict, optional
        Keyword arguments for SmilesGenerator.
    """
    def __init__(self, prefix=None, allow_duplicates=True, **kwargs):
        self.prefix = prefix
        self.allow_duplicates = allow_duplicates
        self.engine = SmilesGenerator(**kwargs)
        self.map = {}

    def add_mol(self, mol):
        """
        Map a molecule name to its corresponding SMILES string and store in the
        SMILES map.

        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        name = mol.GetProp('_Name')
        try:
            int(name)  # check if this is a bare ID
            if self.prefix is None:
                raise TypeError('Bare IDs are not allowed.')
        except ValueError:
            pass
        if self.prefix is not None:
            name = '{}{}'.format(self.prefix, name)
        smiles = self.engine.get_smiles(mol)

        # Failures:
        # * Name is already mapped to a different SMILES
        # * SMILES is already used for a different name
        if name in self.map:  # catch all cases where name is already used
            if self.map[name] != smiles:
                raise ValueError('ID collision for "{}".'.format(name))
        elif not self.allow_duplicates and smiles in self.map.values():
            other = None
            for key, val in self.map.items():
                if val == smiles:
                    other = key
                    break
            raise ValueError(
                'SMILES collision between "{}" and "{}":\n\t{}'.format(
                    name, other, smiles))
        else:
            self.map[name] = smiles

    def get_map(self):
        """
        Get the map.
        """
        return self.map
