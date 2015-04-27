"""
Test NNScore binana featurizer.

TODO(rbharath): Most of these tests are trivial and merely check that the
code can be invoked with trivial arguments. More nontrivial tests require
identification of ligand-receptor structures that boast interesting
geometries (with salt-bridges, pi-cation interactions, etc.)
"""
import unittest

from vs_utils.features.nnscore import NNScoreFeaturizer, binana
from vs_utils.features.nnscore_helper import PDB, atom, point

class TestNNScoreFeaturizer(unittest.TestCase):
  """
  Test NNScoreFeaturizer class.
  """
  def setUp(self):
    """
    Instantiate local featurizer.
    """
    self.featurizer = NNScoreFeaturizer()

class TestBinana(unittest.TestCase):
  """
  Test binana Binding Pose Featurizer.
  """
  def setUp(self):
    """
    Instantiate local copy of binana object.
    """
    self.binana = binana()

    self.receptor = PDB()
    receptor_atom = atom(element="C", coordinates=point(0,0,0))
    self.receptor.AddNewAtom(receptor_atom)

    self.structured_receptor = PDB()
    backbone_alpha_receptor_atom = atom(element="C", coordinates=point(0,0,0),
        structure="ALPHA")
    backbone_beta_receptor_atom = atom(element="C", coordinates=point(0,0,1),
        structure="BETA")
    backbone_other_receptor_atom = atom(element="C", coordinates=point(0,0,2),
        structure="OTHER")
    sidechain_alpha_receptor_atom = atom(element="S", coordinates=point(0,0,3),
        structure="ALPHA")
    sidechain_beta_receptor_atom = atom(element="S", coordinates=point(0,0,4),
        structure="BETA")
    sidechain_other_receptor_atom = atom(element="S", coordinates=point(0,0,5),
        structure="OTHER")
    self.structured_receptor.AddNewAtoms([backbone_alpha_receptor_atom,
        backbone_beta_receptor_atom, backbone_other_receptor_atom,
        sidechain_alpha_receptor_atom, sidechain_beta_receptor_atom,
        sidechain_other_receptor_atom])
    

    self.hydrophobic_ligand = PDB()
    hydrophobic_ligand_atom = atom(element="C", coordinates=point(0,0,1))
    self.hydrophobic_ligand.AddNewAtom(hydrophobic_ligand_atom)

  def testComputeHydrophobicContact(self):
    """
    TestBinana: Test that hydrophobic contacts are established.
    """
    hydrophobics = self.binana.compute_hydrophobic_contacts(
        self.hydrophobic_ligand, self.receptor)
    #assert len(hydrophobics.keys()) == 1

  def testComputeElectrostaticEnergy(self):
    """
    TestBinana: Test that electrostatic energies are computed.
    """
    ligand_receptor_electrostatics = (
        self.binana.compute_electrostatic_energy(
            self.hydrophobic_ligand, self.receptor))
    # TODO(bramsundar): Add a more nontrivial test of electrostatics here.
    #assert len(ligand_receptor_electrostatics) == 1

  def testComputeActiveSiteFlexibility(self):
    """
    TestBinana: Gather statistics about active site protein atoms.
    """
    active_site_flexibility = (
        self.binana.compute_active_site_flexibility(self.hydrophobic_ligand,
            self.structured_receptor))
    print active_site_flexibility
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
      self.binana.compute_hydrogen_bonds(self.hydrophobic_ligand,
          self.receptor))

  def testComputeLigandAtomCounts(self):
    """
    TestBinana: Compute the Number of Ligand Atom Counts.
    """
    ligand_atom_types = (
      self.binana.compute_ligand_atom_counts(self.hydrophobic_ligand))

  def testComputeLigandReceptorContacts(self):
    """
    TestBinana: Compute contacts between Ligand and receptor.
    """
    ligand_receptor_close_contacts, ligand_receptor_contacts = (
      self.binana.compute_ligand_receptor_contacts(self.hydrophobic_ligand,
          self.receptor))

  def testComputePiPiStacking(self):
    """
    TestBinana: Compute Pi-Pi Stacking.
    """
    # TODO(bramsundar): THERE ARE NO AROMATIC RINGS HERE! BOGUS TEST!
    pi_stacking = (
        self.binana.compute_pi_pi_stacking(self.hydrophobic_ligand,
            self.receptor))


  def testComputePiT(self):
    """
    TestBinana: Compute Pi-T Interactions.
    """
    # TODO(bramsundar): THERE ARE NO AROMATIC RINGS HERE! BOGUS TEST!
    pi_T = (
        self.binana.compute_pi_T(self.hydrophobic_ligand,
            self.receptor))

  def testComputePiCation(self):
    """
    TestBinana: Compute Pi-Cation Interactions.
    """
    # TODO(bramsundar): THERE ARE NO AROMATIC RINGS HERE! BOGUS TEST!
    pi_cation = (
        self.binana.compute_pi_cation(self.hydrophobic_ligand,
            self.receptor))

  def testComputeSaltBridges(self):
    """
    TestBinana: Compute Salt Bridges.
    """
    # TODO(bramsundar): THERE ARE NO AROMATIC RINGS HERE! BOGUS TEST!
    salt_bridges = (
        self.binana.compute_salt_bridges(self.hydrophobic_ligand,
            self.receptor))

  def testComputeInputVector(self):
    """
    TestBinana: Compute Input Vector.
    """
    input_vector = (
        self.binana.compute_input_vector(self.hydrophobic_ligand,
            self.receptor))
