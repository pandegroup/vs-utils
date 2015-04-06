"""
Test NNScore binana featurizer.
"""
import unittest

from vs_utils.features.nnscore import NNScoreFeaturizer, binana, command_line_parameters
from vs_utils.features.nnscore_helper import PDB, atom

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
    receptor_atom = atom(element="C")
    self.receptor.AddNewAtom(receptor_atom)

    self.hydrophobic_ligand = PDB()
    hydrophobic_ligand_atom = atom(element="C")
    self.hydrophobic_ligand.AddNewAtom(hydrophobic_ligand_atom)

  def testComputeHydrophobicContact(self):
    """
    TestBinana: Test that hydrophobic contacts are established.
    """
    hydrophobics, pdb_hydrophobic = self.binana.compute_hydrophobic_contact(
        self.hydrophobic_ligand, self.receptor)
    assert len(pdb_hydrophobic.AllAtoms.keys()) == 2
    assert len(hydrophobics.keys()) == 1
