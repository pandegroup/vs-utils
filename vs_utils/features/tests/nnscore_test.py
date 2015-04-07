"""
Test NNScore binana featurizer.
"""
import unittest

from vs_utils.features.nnscore import NNScoreFeaturizer, binana, command_line_parameters
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
    backbone_beta_receptor_atom = atom(element="C", coordinates=point(0,0,0),
        structure="BETA")
    backbone_other_receptor_atom = atom(element="C", coordinates=point(0,0,0),
        structure="OTHER")
    sidechain_alpha_receptor_atom = atom(element="S", coordinates=point(0,0,0),
        structure="ALPHA")
    sidechain_beta_receptor_atom = atom(element="S", coordinates=point(0,0,0),
        structure="BETA")
    sidechain_other_receptor_atom = atom(element="S", coordinates=point(0,0,0),
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
    hydrophobics, pdb_hydrophobic = self.binana.compute_hydrophobic_contact(
        self.hydrophobic_ligand, self.receptor)
    assert len(pdb_hydrophobic.AllAtoms.keys()) == 2
    assert len(hydrophobics.keys()) == 1

  def testComputeElectrostaticEnergy(self):
    """
    TestBinana: Test that electrostatic energies are computed.
    """
    ligand_receptor_electrostatics = (
        self.binana.compute_electrostatic_energy(
            self.hydrophobic_ligand, self.receptor))
    # TODO(bramsundar): Add a more nontrivial test of electrostatics here.
    assert len(ligand_receptor_electrostatics) == 1

  def testComputeActiveSiteFlexibility(self):
    """
    TestBinana: Gather statistics about active site protein atoms.
    """
    active_size_flexibility = (
        self.binana.compute_active_site_flexibility(self.hydrophobic_ligand,
        self.structured_receptor))
