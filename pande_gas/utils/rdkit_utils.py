"""
RDKit utility functions.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import gzip
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem


def read(filename, mol_format=None, return_generator=False):
    """
    Read molecules. Supports SDF or SMILES format.

    Parameters
    ----------
    filename : str
        Filename.
    mol_format : str, optional
        Molecule file format. Currently supports 'sdf' and 'smi'.
    return_generator : bool, optional (default False)
        Whether to return a Mol generator. If False, all molecules are read
        into memory and an ndarray is returned.
    """
    mols = read_generator(filename, mol_format)
    if return_generator:
        return mols
    else:
        mols = np.asarray(list(mols))
        return mols


def read_generator(filename, mol_format=None):
    """
    Read molecules. Supports SDF or SMILES format.

    Parameters
    ----------
    filename : str
        Filename.
    mol_format : str, optional
        Molecule file format. Currently supports 'sdf' and 'smi'.
    """
    if mol_format is None:
        if filename.endswith(('.sdf', '.sdf.gz')):
            mol_format = 'sdf'
        elif filename.endswith(('.smi', '.smi.gz', '.can', '.can.gz', '.ism',
                                '.ism.gz')):
            mol_format = 'smi'
        else:
            raise NotImplementedError('Unable to guess molecule file format.')
    if mol_format not in ['sdf', 'smi']:
        raise NotImplementedError('Unsupported molecule file format ' +
                                  '"{}".'.format(mol_format))
    if filename.endswith('.gz'):
        f = gzip.open(filename)
    else:
        f = open(filename)
    if mol_format == 'sdf':
        for mol in Chem.ForwardSDMolSupplier(f):
            yield mol
    elif mol_format == 'smi':
        lines = [line.strip().split() for line in f.readlines()]
        for line in lines:
            if len(line) > 1:
                smiles, name = line
            else:
                smiles = line
                name = None
            mol = Chem.MolFromSmiles(smiles)
            if name is not None:
                mol.SetProp('_Name', name)
            yield mol
    f.close()


def write(mols, filename):
    """
    Write SDF molecules.

    Parameters
    ----------
    mols : list
        Molecules to write.
    filename : str
        Output filename.
    """
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'wb')
    else:
        f = open(filename, 'wb')
    w = Chem.SDWriter(f)
    for mol in np.atleast_1d(mols):
        if mol.GetNumConformers():
            for conf in mol.GetConformers():
                w.write(mol, confId=conf.GetId())
        else:
            w.write(mol)
    w.close()
    f.close()


def generate_conformers(mol, n_conformers=1, rmsd_threshold=0.5):
    """
    Generate molecule conformers. See:
    * http://rdkit.org/docs/GettingStartedInPython.html
      #working-with-3d-molecules
    * http://pubs.acs.org/doi/full/10.1021/ci2004658

    This function returns a copy of the molecule, created before adding
    hydrogens.

    Parameters
    ----------
    mol : RDKit Mol
        Molecule.
    n_conformers : int, optional (default 1)
        Maximum number of conformers to generate.
    rmsd_threshold : float, optional (default 0.5)
        RMSD threshold for distinguishing conformers.
    """
    mol = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol, n_conformers,
                                      pruneRmsThresh=rmsd_threshold)
    assert mol.GetNumConformers() >= 1
    cids = np.asarray(cids, dtype=int)

    # minimize conformers and get energies
    energy = np.zeros(cids.size, dtype=float)
    for cid in cids:
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
        ff.Minimize()
        energy[cid] = ff.CalcEnergy()

    # calculate RMSD between minimized conformers
    rmsd = np.zeros((cids.size, cids.size), dtype=float)
    for i in xrange(cids.size):
        for j in xrange(cids.size):
            if i >= j:
                continue
            rmsd[i, j] = AllChem.GetBestRMS(mol, mol, cids[i], cids[j])
            rmsd[j, i] = rmsd[i, j]

    # discard conformers within RMSD threshold
    _, discard = choose_conformers(energy, rmsd, rmsd_threshold)
    for i in discard:
        mol.RemoveConformer(cids[i])
    return mol


def choose_conformers(energy, rmsd, rmsd_threshold=0.5):
    """
    Select diverse conformers starting with lowest energy.

    Parameters
    ----------
    energy : list
        Conformer energies.
    rmsd : ndarray
        Conformer-conformer RMSD values.
    rmsd_threshold : float, optional (default 0.5)
        RMSD threshold for distinguishing conformers.
    """
    if len(energy) == 1:
        return [0], []
    sort = np.argsort(energy)
    keep = [sort[0]]
    discard = []
    for i in sort[1:]:
        this_rmsd = rmsd[i][np.asarray(keep, dtype=int)]
        if np.all(this_rmsd >= rmsd_threshold):
            keep.append(i)
        else:
            discard.append(i)
    keep = np.asarray(keep, dtype=int)
    discard = np.asarray(discard, dtype=int)
    return keep, discard


def interatomic_distances(conf):
    """
    Get interatomic distances for atoms in a molecular conformer.

    Parameters
    ----------
    conf : RDKit Conformer
        Molecule conformer.
    """
    n_atoms = conf.GetNumAtoms()
    coords = [conf.GetAtomPosition(i) for i in xrange(n_atoms)]
    d = np.zeros((n_atoms, n_atoms), dtype=float)
    for i in xrange(n_atoms):
        for j in xrange(n_atoms):
            if i < j:
                d[i, j] = coords[i].Distance(coords[j])
                d[j, i] = d[i, j]
            else:
                continue
    return d
