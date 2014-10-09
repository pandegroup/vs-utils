"""
Dataset preparation utilities.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import gzip

from rdkit import Chem

from . import SmilesGenerator


class MoleculeDatabase(object):
    """
    Molecule database.

    Molecules are keyed by SMILES.

    Parameters
    ----------
    kwargs : dict, optional
        Keyword arguments for SmilesMap.
    """
    def __init__(self, **kwargs):
        self.engine = SmilesGenerator(**kwargs)
        self.smiles = set()

    def __len__(self):
        return len(self.smiles)

    def __iter__(self):
        for smiles in self.smiles:
            yield smiles

    def add_mol(self, mol):
        """
        Add a molecule to the database.

        Parameters
        ----------
        mol : RDKit Mol
            Molecule.
        """
        self.smiles.add(self.engine.get_smiles(mol))

    def load(self, filename):
        """
        Load an existing database.

        Parameters
        ----------
        filename : str
            Existing database filename.
        """
        if filename.endswith('.gz'):
            f = gzip.open(filename)
        else:
            f = open(filename)
        for line in f:
            smiles = line.strip()
            mol = Chem.MolFromSmiles(smiles)  # sanity check
            if mol is None:
                raise ValueError(
                    'Database is unreadable: "{}".'.format(smiles))
            self.smiles.add(smiles)
        f.close()

    def save(self, filename):
        """
        Save the database to disk.

        Parameters
        ----------
        filename : str
            Filename.
        """
        if filename.endswith('.gz'):
            f = gzip.open(filename, 'wb')
        else:
            f = open(filename, 'wb')
        for smiles in self.smiles:
            f.write('{}\n'.format(smiles))
        f.close()
