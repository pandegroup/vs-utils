"""
RDKit utility functions.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

from collections import OrderedDict
import gzip
import numpy as np
import warnings

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.SaltRemover import SaltRemover


def read_mols(filename, mol_format=None, remove_salts=True):
    """
    Read molecules. Supports SDF or SMILES format.

    Molecular conformers are combined into a single molecule. The criteria
    for conformers are:
    * Identical SMILES strings
    * Identical compound names
    * Contiguous in the file

    Parameters
    ----------
    filename : str
        Filename.
    mol_format : str, optional
        Molecule file format. Currently supports 'sdf' and 'smi'.
    remove_salts : bool, optional (default True)
        Whether to remove salts.

    Returns
    -------
    An ndarray containing RDMol objects.
    """
    mols = []
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
            mols.append(mol)
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
            mols.append(mol)
    f.close()

    # remove salts
    if remove_salts:
        salt_remover = SaltRemover()
        mols = [salt_remover.StripMol(mol) for mol in mols]

    # combine conformers
    smiles = np.asarray([Chem.MolToSmiles(mol, isomericSmiles=True,
                                          canonical=True) for mol in mols])
    names = []
    for mol in mols:
        if mol.HasProp('_Name'):
            names.append(mol.GetProp('_Name'))
        else:
            names.append(None)
    to_combine = [[mols[0]]]
    idx = 0
    for i in xrange(1, smiles.size):
        if smiles[i] == smiles[i - 1] and names[i] == names[i - 1]:
            if names[i] is None:
                warnings.warn("Combining conformers for unnamed molecules.")
            to_combine[idx].append(mols[i])
        else:
            to_combine.append([mols[i]])
            idx += 1

    combined = []
    for mols in to_combine:
        if len(mols) == 1:
            mol, = mols
        else:
            mol = Chem.Mol(mols[0])
            for other in mols[1:]:
                for conf in other.GetConformers():
                    mol.AddConformer(conf)
        combined.append(mol)
    combined = np.asarray(combined)
    return combined


def write_mols(mols, filename):
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
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=int(cid))
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
