"""
Miscellaneous utility functions.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import numpy as np
import os

from rdkit import Chem
from rdkit_utils import PicklableMol, serial


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


class SmilesMap(object):
    """
    Map compound names to SMILES.

    Parameters
    ----------
    prefix : str, optional
        Prefix to prepend to IDs.
    """
    def __init__(self, prefix=None):
        self.prefix = prefix
        self.map = {}

    def add_mol(self, mol):
        """
        Map a molecule name to its corresponding SMILES string.

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
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        if name in self.map and self.map[name] != smiles:
            raise ValueError('ID collision for "{}".'.format(name))
        else:
            self.map[name] = smiles

    def get_map(self):
        """
        Get the map.
        """
        return self.map
