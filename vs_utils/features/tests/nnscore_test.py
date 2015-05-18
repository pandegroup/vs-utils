"""
Test NNScore Binana featurizer.

TODO(rbharath): Most of these tests are simplistic. More nontrivial tests
require identification of ligand-receptor structures that boast interesting
interactions (with salt-bridges, pi-cation interactions, etc.)
"""
import os
import numpy as np
import re
import unittest
import itertools

from vs_utils.features.nnscore import Binana
from vs_utils.utils.nnscore_pdb import PDB
from vs_utils.utils.nnscore_utils import Atom
from vs_utils.utils.nnscore_utils import Point
from vs_utils.utils.tests import __file__ as test_directory

def data_dir():
  """Get location of data directory."""
  return os.path.join(os.path.dirname(test_directory), "data")

class TestBinana(unittest.TestCase):
  """
  Test Binana Binding Pose Featurizer.
  """
  def setUp(self):
    """
    Instantiate local copy of Binana object.
    """
    self.binana = Binana()

    self.prgr_receptor = PDB()
    prgr_pdb = os.path.join(data_dir(), "prgr_hyd.pdb")
    prgr_pdbqt = os.path.join(data_dir(), "prgr_hyd.pdbqt")
    self.prgr_receptor.load_from_files(prgr_pdb, prgr_pdbqt)

    # This compound is CHEMBL1164248
    self.prgr_active = PDB()
    prgr_active_pdb = os.path.join(data_dir(), "prgr_active0_hyd.pdb")
    prgr_active_pdbqt = os.path.join(data_dir(), "prgr_active0_hyd.pdbqt")
    self.prgr_active.load_from_files(prgr_active_pdb, prgr_active_pdbqt)

    self.cAbl_receptor = PDB()
    cAbl_receptor_pdb = os.path.join(data_dir(), "c-Abl_hyd.pdb")
    cAbl_receptor_pdbqt = os.path.join(data_dir(), "c-Abl_hyd.pdbqt")
    self.cAbl_receptor.load_from_files(cAbl_receptor_pdb,
        cAbl_receptor_pdbqt)

    # This compound is imatinib
    self.cAbl_active = PDB()
    cAbl_active_pdb = os.path.join(data_dir(), "imatinib_hyd.pdb")
    cAbl_active_pdbqt = os.path.join(data_dir(), "imatinib_hyd.pdbqt")
    self.cAbl_active.load_from_files(cAbl_active_pdb,
        cAbl_active_pdbqt)

  def testComputeHydrophobicContacts(self):
    """
    TestBinana: Test that hydrophobic contacts are established.
    """
    prgr_hydrophobics = self.binana.compute_hydrophobic_contacts(
        self.prgr_active, self.prgr_receptor)
    cAbl_hydrophobics = self.binana.compute_hydrophobic_contacts(
        self.cAbl_active, self.cAbl_receptor)
    for hydrophobics in [prgr_hydrophobics, cAbl_hydrophobics]:
      assert len(hydrophobics) == 6
      assert "BACKBONE_ALPHA" in hydrophobics
      assert "BACKBONE_BETA" in hydrophobics
      assert "BACKBONE_OTHER" in hydrophobics
      assert "SIDECHAIN_ALPHA" in hydrophobics
      assert "SIDECHAIN_BETA" in hydrophobics
      assert "SIDECHAIN_OTHER" in hydrophobics

  def testComputeElectrostaticEnergy(self):
    """
    TestBinana: Test that electrostatic energies are computed.

    TODO(rbharath): This method reports that no electrostatics are found
    for either prgr or cAbl. I'm pretty sure this is a bug.
    """
    prgr_electro = (
        self.binana.compute_electrostatic_energy(
            self.prgr_active, self.prgr_receptor))
    cAbl_electro = (
        self.binana.compute_electrostatic_energy(
            self.cAbl_active, self.cAbl_receptor))
    for electrostatics in [prgr_electro, cAbl_electro]:
      # The keys of these dicts are pairs of atomtypes, but the keys are
      # sorted so that ("C", "O") is always written as "C_O". Thus, for N
      # atom types, there are N*(N+1)/2 unique pairs.
      N = len(Binana.atom_types)
      assert len(electrostatics) == N*(N+1)/2
      assert np.count_nonzero(np.array(electrostatics.values())) > 0

  def testComputeActiveSiteFlexibility(self):
    """
    TestBinana: Gather statistics about active site protein atoms.
    """
    prgr_flex = (
        self.binana.compute_active_site_flexibility(self.prgr_active,
            self.prgr_receptor))
    cAbl_flex = (
        self.binana.compute_active_site_flexibility(self.cAbl_active,
            self.cAbl_receptor))
    for active_site_flexibility in [prgr_flex, cAbl_flex]:
      assert len(active_site_flexibility.keys()) == 6
      assert "BACKBONE_ALPHA" in active_site_flexibility
      assert "BACKBONE_BETA" in active_site_flexibility
      assert "BACKBONE_OTHER" in active_site_flexibility
      assert "SIDECHAIN_ALPHA" in active_site_flexibility
      assert "SIDECHAIN_BETA" in active_site_flexibility
      assert "SIDECHAIN_OTHER" in active_site_flexibility

  def testComputeHydrogenBonds(self):
    """
    TestBinana: Compute the number of hydrogen bonds.

    TODO(rbharath): The hydrogen-bond angle cutoff seems like it's
    incorrect to me. The hydrogens are placed by openbabel and aren't
    optimized, so I'm pretty sure that this code will miss many hydrogens.
    Here are some options:
    -) Find a method to optimize the hydrogen placement.
    -) Place a more permissive angle cutoff for hydrogens.
    -) Allow for "buckets": angles 0-20, 20-40, 40-60, etc. and count the
    number of hydrogen bonds in each bucket.
    """
    prgr_hbonds = (
      self.binana.compute_hydrogen_bonds(self.prgr_active,
          self.prgr_receptor))
    cAbl_hbonds = (
      self.binana.compute_hydrogen_bonds(self.cAbl_active,
          self.cAbl_receptor))
    for hbonds in [prgr_hbonds, cAbl_hbonds]:
      assert len(hbonds) == 12
      assert "HDONOR-LIGAND_BACKBONE_ALPHA" in hbonds
      assert "HDONOR-LIGAND_BACKBONE_BETA" in hbonds
      assert "HDONOR-LIGAND_BACKBONE_OTHER" in hbonds
      assert "HDONOR-LIGAND_SIDECHAIN_ALPHA" in hbonds
      assert "HDONOR-LIGAND_SIDECHAIN_BETA" in hbonds
      assert "HDONOR-LIGAND_SIDECHAIN_OTHER" in hbonds
      assert "HDONOR-RECEPTOR_BACKBONE_ALPHA" in hbonds
      assert "HDONOR-RECEPTOR_BACKBONE_BETA" in hbonds
      assert "HDONOR-RECEPTOR_BACKBONE_OTHER" in hbonds
      assert "HDONOR-RECEPTOR_SIDECHAIN_ALPHA" in hbonds
      assert "HDONOR-RECEPTOR_SIDECHAIN_BETA" in hbonds
      assert "HDONOR-RECEPTOR_SIDECHAIN_OTHER" in hbonds

  def testComputeLigandAtomCounts(self):
    """
    TestBinana: Compute the Number of Ligand Atom Counts.
    """
    prgr_counts = (
      self.binana.compute_ligand_atom_counts(self.prgr_active))
    cAbl_counts = (
      self.binana.compute_ligand_atom_counts(self.cAbl_active))
    for ligand_atom_counts in [prgr_counts, cAbl_counts]:
      assert len(ligand_atom_counts) == len(Binana.atom_types)

  def testComputeLigandReceptorContacts(self):
    """
    TestBinana: Compute contacts between Ligand and receptor.
    """
    prgr_close, prgr_contacts = (
      self.binana.compute_ligand_receptor_contacts(self.prgr_active,
          self.prgr_receptor))
    cAbl_close, cAbl_contacts = (
      self.binana.compute_ligand_receptor_contacts(self.cAbl_active,
          self.cAbl_receptor))
    for ligand_receptor_close_contacts, ligand_receptor_contacts in [
        (prgr_close, prgr_contacts), (cAbl_close, cAbl_contacts)]:
      print "close_contacts"
      for key in ligand_receptor_close_contacts:
        val = ligand_receptor_close_contacts[key]
        if val != 0:
          print (key, val)
      print "contacts"
      for key in ligand_receptor_contacts:
        val = ligand_receptor_contacts[key]
        if val != 0:
          print (key, val)
      # The keys of these dicts are pairs of atomtypes, but the keys are
      # sorted so that ("C", "O") is always written as "C_O". Thus, for N
      # atom types, there are N*(N+1)/2 unique pairs.
      N = len(Binana.atom_types)
      assert len(ligand_receptor_close_contacts) == N*(N+1)/2
      assert len(ligand_receptor_contacts) == N*(N+1)/2

  def testComputePiPiStacking(self):
    """
    TestBinana: Compute Pi-Pi Stacking.
    """
    # TODO(bramsundar): prgr has no pi-pi stacking. Find a different
    # complex that does.
    prgr_pi_stacking = (
        self.binana.compute_pi_pi_stacking(self.prgr_active,
            self.prgr_receptor))
    cAbl_pi_stacking = (
        self.binana.compute_pi_pi_stacking(self.cAbl_active,
            self.cAbl_receptor))
    for pi_stacking in [prgr_pi_stacking, cAbl_pi_stacking]:
      assert len(pi_stacking) == 3
      assert "STACKING_ALPHA" in pi_stacking
      assert "STACKING_BETA" in pi_stacking
      assert "STACKING_OTHER" in pi_stacking


  def testComputePiT(self):
    """
    TestBinana: Compute Pi-T Interactions.

    TODO(rbharath): I believe that the imatininb-cAbl complex has a pi-T
    interaction. This code has a bug since it reports that no such
    interaction is found.
    """
    # TODO(bramsundar): prgr has no pi-T interactions. Find an alternative
    # structure that does.
    prgr_pi_T = (
        self.binana.compute_pi_T(self.prgr_active,
            self.prgr_receptor))
    cAbl_pi_T = (
        self.binana.compute_pi_T(self.cAbl_active,
            self.cAbl_receptor))
    for pi_T in [prgr_pi_T, cAbl_pi_T]:
      assert len(pi_T) == 3
      assert "T-SHAPED_ALPHA" in pi_T
      assert "T-SHAPED_BETA" in pi_T
      assert "T-SHAPED_OTHER" in pi_T

  def testComputePiCation(self):
    """
    TestBinana: Compute Pi-Cation Interactions.
    """
    # TODO(rbharath): prgr doesn't have any pi-cation interactions. Find a
    # different complex that exhibits this interaction.
    prgr_pi_cation = (
        self.binana.compute_pi_cation(self.prgr_active,
            self.prgr_receptor))
    cAbl_pi_cation = (
        self.binana.compute_pi_cation(self.cAbl_active,
            self.cAbl_receptor))
    for pi_cation in [prgr_pi_cation, cAbl_pi_cation]:
      assert len(pi_cation) == 6
      assert 'PI-CATION_LIGAND-CHARGED_ALPHA' in pi_cation
      assert 'PI-CATION_LIGAND-CHARGED_BETA' in pi_cation
      assert 'PI-CATION_LIGAND-CHARGED_OTHER' in pi_cation
      assert 'PI-CATION_RECEPTOR-CHARGED_ALPHA' in pi_cation
      assert 'PI-CATION_RECEPTOR-CHARGED_BETA' in pi_cation
      assert 'PI-CATION_RECEPTOR-CHARGED_OTHER' in pi_cation

  def testComputeSaltBridges(self):
    """
    TestBinana: Compute Salt Bridges.
    """
    # TODO(bramsundar): prgr contains no salt-bridge interactions. Find a
    # complex with an actual salt-bridge interaction.
    prgr_bridges = (
        self.binana.compute_salt_bridges(self.prgr_active,
            self.prgr_receptor))
    cAbl_bridges = (
        self.binana.compute_salt_bridges(self.cAbl_active,
            self.cAbl_receptor))
    for salt_bridges in [prgr_bridges, cAbl_bridges]:
      assert len(salt_bridges) == 3
      assert 'SALT-BRIDGE_ALPHA' in salt_bridges
      assert 'SALT-BRIDGE_BETA' in salt_bridges
      assert 'SALT-BRIDGE_OTHER' in salt_bridges


  def testComputeInputVector(self):
    """
    TestBinana: Compute Input Vector.
    """
    prgr_vector = (
        self.binana.compute_input_vector(self.prgr_active,
            self.prgr_receptor))
    cAbl_vector = (
        self.binana.compute_input_vector(self.cAbl_active,
            self.cAbl_receptor))
    N = len(Binana.atom_types)
    # Lengths:
    # ligand_receptor_close_contacts: N*(N+1)/2
    # ligand_receptor_contacts: N*(N+1)/2
    # ligand_receptor_electrostatics: N*(N+1)/2
    # ligand_atom_counts: N
    # hbonds: 12
    # hydrophobics: 6
    # stacking: 3
    # pi_cation: 6
    # t_shaped: 3
    # active_site_flexibility: 6
    # salt_bridges: 3
    # rotatable_boonds_count: 1
    total_len = 3*N*(N+1)/2 + N + 12 + 6 + 3 + 6 + 3 + 6 + 3 + 1
    for input_vector in [prgr_vector, cAbl_vector]:
      assert len(input_vector) == total_len
