"""
Test NNScore binana featurizer.
"""
import unittest

from vs_utils.features.nnscore import NNScoreFeaturizer, binana, command_line_parameters

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
    self.parameters = command_line_parameters()
