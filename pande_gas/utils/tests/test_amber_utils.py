"""
Tests for amber_utils.
"""
from cStringIO import StringIO
import numpy as np

from rdkit import Chem
from rdkit_utils import conformers

from pande_gas.utils import amber_utils, pdb_utils


def test_antechamber_charges_and_radii():
    """Test Antechamber charges and radii."""
    mol = Chem.MolFromSmiles(test_smiles)
    mol = conformers.generate_conformers(mol)
    antechamber = amber_utils.Antechamber()
    charges, radii = antechamber.get_charges_and_radii(mol)
    assert charges.size == 12  # 12 atoms, 6 carbon, 6 hydrogen
    assert np.allclose(charges.sum(), 0)  # neutral molecule
    assert radii.size == 12  # 12 atoms
    assert np.count_nonzero(radii > 0)  # no zero radii


def test_pbsa_esp_grid():
    """Test PBSA ESP grid."""
    mol = Chem.MolFromSmiles(test_smiles)
    mol = conformers.generate_conformers(mol)

    # generate PQR
    reader = pdb_utils.PdbReader()
    antechamber = amber_utils.Antechamber()
    charges, radii = antechamber.get_charges_and_radii(mol)
    pdb = Chem.MolToPDBBlock(mol)
    pqr = reader.pdb_to_pqr(StringIO(pdb), charges, radii)

    # calculate ESP grid
    pbsa = amber_utils.PBSA()
    grid, _ = pbsa.get_esp_grid(pqr)

    # the grid should be cubic
    size = grid.shape[0]
    assert grid.shape == (size, size, size), grid.shape

    # and not be all zeros
    assert np.count_nonzero(grid)

test_smiles = 'c1ccccc1'
