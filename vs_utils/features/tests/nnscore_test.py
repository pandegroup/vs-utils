"""
Test NNScore binana featurizer.

TODO(rbharath): Most of these tests are trivial and merely check that the
code can be invoked with trivial arguments. More nontrivial tests require
identification of ligand-receptor structures that boast interesting
geometries (with salt-bridges, pi-cation interactions, etc.)
"""
import os
import unittest
import itertools

from vs_utils.features.nnscore import binana
from vs_utils.features.nnscore_helper import PDB, Atom, Point
from vs_utils.features.tests import __file__ as test_directory

class TestNNScoreFeaturizer(unittest.TestCase):
  """
  Test NNScoreFeaturizer class.
  """
  def setUp(self):
    """
    Instantiate local featurizer.
    """
    self.featurizer = NNScoreFeaturizer()

def data_dir():
  """Get location of data directory."""
  return os.path.join(os.path.dirname(test_directory), "data")

class TestBinana(unittest.TestCase):
  """
  Test binana Binding Pose Featurizer.
  """
  def setUp(self):
    """
    Instantiate local copy of binana object.
    """
    self.binana = binana()

    self.prgr_receptor = PDB()
    prgr_receptor_path = os.path.join(data_dir(), "prgr.pdb")
    self.prgr_receptor.load_PDB_from_file(prgr_receptor_path)

    # TODO(rbharath): Document the pubchem id of this active.
    self.prgr_active = PDB()
    prgr_active_path = os.path.join(data_dir(), "prgr_active0.pdb")
    self.prgr_active.load_PDB_from_file(prgr_active_path)

  def testComputeHydrophobicContact(self):
    """
    TestBinana: Test that hydrophobic contacts are established.
    """
    hydrophobics = self.binana.compute_hydrophobic_contacts(
        self.prgr_active, self.prgr_receptor)
    #assert len(hydrophobics.keys()) == 1

  def testComputeElectrostaticEnergy(self):
    """
    TestBinana: Test that electrostatic energies are computed.
    """
    ligand_receptor_electrostatics = (
        self.binana.compute_electrostatic_energy(
            self.prgr_active, self.prgr_receptor))
    # TODO(bramsundar): Add a more nontrivial test of electrostatics here.
    # The keys of these dicts are pairs of atomtypes, but the keys are
    # sorted so that ("C", "O") is always written as "C_O". Thus, for N
    # atom types, there are N*(N+1)/2 unique pairs.
    N = len(binana.atom_types)
    print N
    print ligand_receptor_electrostatics
    print len(ligand_receptor_electrostatics)
    print N*(N+1)/2
    for first, second in itertools.product(binana.atom_types,
      binana.atom_types):
      key = "_".join(sorted([first, second]))
      ligand_receptor_electrostatics.pop(key, None)
    print ligand_receptor_electrostatics
    assert len(ligand_receptor_electrostatics) == N*(N+1)/2

  def testComputeActiveSiteFlexibility(self):
    """
    TestBinana: Gather statistics about active site protein atoms.
    """
    active_site_flexibility = (
        self.binana.compute_active_site_flexibility(self.prgr_active,
            self.prgr_receptor))
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
    """
    # TODO(bramsundar): Add a nontrivial test here
    hbonds = (
      self.binana.compute_hydrogen_bonds(self.prgr_active,
          self.prgr_receptor))
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
    ligand_atom_counts = (
      self.binana.compute_ligand_atom_counts(self.prgr_active))
    assert len(ligand_atom_counts) == len(binana.atom_types)
    assert ligand_atom_counts["A"] == 0
    assert ligand_atom_counts["BR"] == 0
    assert ligand_atom_counts["C"] == 27
    assert ligand_atom_counts["CD"] == 0
    assert ligand_atom_counts["CL"] == 0
    assert ligand_atom_counts["CU"] == 0
    assert ligand_atom_counts["F"] == 0
    assert ligand_atom_counts["FE"] == 0
    assert ligand_atom_counts["H"] == 1
    assert ligand_atom_counts["HD"] == 0
    assert ligand_atom_counts["I"] == 0
    assert ligand_atom_counts["MG"] == 0
    assert ligand_atom_counts["MN"] == 0
    assert ligand_atom_counts["N"] == 2
    assert ligand_atom_counts["NA"] == 0
    assert ligand_atom_counts["O"] == 3
    assert ligand_atom_counts["OA"] == 0
    assert ligand_atom_counts["P"] == 0
    assert ligand_atom_counts["S"] == 0
    assert ligand_atom_counts["SA"] == 0
    assert ligand_atom_counts["ZN"] == 0

  def testComputeLigandReceptorContacts(self):
    """
    TestBinana: Compute contacts between Ligand and receptor.
    """
    ligand_receptor_close_contacts, ligand_receptor_contacts = (
      self.binana.compute_ligand_receptor_contacts(self.prgr_active,
          self.prgr_receptor))
    # The keys of these dicts are pairs of atomtypes, but the keys are
    # sorted so that ("C", "O") is always written as "C_O". Thus, for N
    # atom types, there are N*(N+1)/2 unique pairs.
    N = len(binana.atom_types)
    print N
    print ligand_receptor_close_contacts
    print len(ligand_receptor_close_contacts)
    print N*(N+1)/2
    for first, second in itertools.product(binana.atom_types,
      binana.atom_types):
      key = "_".join(sorted([first, second]))
      ligand_receptor_close_contacts.pop(key, None)
    print ligand_receptor_close_contacts
    #assert len(ligand_receptor_close_contacts) == N*(N+1)/2
    assert 0 == 1

  def testComputePiPiStacking(self):
    """
    TestBinana: Compute Pi-Pi Stacking.
    """
    # TODO(bramsundar): prgr has no pi-pi stacking. Find a different
    # complex that does. 
    pi_stacking = (
        self.binana.compute_pi_pi_stacking(self.prgr_active,
            self.prgr_receptor))
    assert len(pi_stacking) == 3
    assert "STACKING_ALPHA" in pi_stacking
    assert "STACKING_BETA" in pi_stacking
    assert "STACKING_OTHER" in pi_stacking


  def testComputePiT(self):
    """
    TestBinana: Compute Pi-T Interactions.
    """
    # TODO(bramsundar): prgr has no pi-T interactions. Find an alternative
    # structure that does. 
    pi_T = (
        self.binana.compute_pi_T(self.prgr_active,
            self.prgr_receptor))
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
    pi_cation = (
        self.binana.compute_pi_cation(self.prgr_active,
            self.prgr_receptor))
    assert len(pi_cation) == 6
    assert 'LIGAND-CHARGED_ALPHA' in pi_cation
    assert 'LIGAND-CHARGED_BETA' in pi_cation
    assert 'LIGAND-CHARGED_OTHER' in pi_cation
    assert 'RECEPTOR-CHARGED_ALPHA' in pi_cation
    assert 'RECEPTOR-CHARGED_BETA' in pi_cation
    assert 'RECEPTOR-CHARGED_OTHER' in pi_cation

  def testComputeSaltBridges(self):
    """
    TestBinana: Compute Salt Bridges.
    """
    # TODO(bramsundar): prgr contains no salt-bridge interactions. Find a
    # complex with an actual salt-bridge interaction.
    salt_bridges = (
        self.binana.compute_salt_bridges(self.prgr_active,
            self.prgr_receptor))
    assert len(salt_bridges) == 3
    assert 'SALT-BRIDGE_ALPHA' in salt_bridges
    assert 'SALT-BRIDGE_BETA' in salt_bridges
    assert 'SALT-BRIDGE_OTHER' in salt_bridges


  def testComputeInputVector(self):
    """
    TestBinana: Compute Input Vector.
    """
    input_vector = (
        self.binana.compute_input_vector(self.prgr_active,
            self.prgr_receptor))
