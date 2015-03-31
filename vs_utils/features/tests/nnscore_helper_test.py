"""
Test NNScore Helper Classes
"""
import unittest
import numpy as np

from vs_utils.features.nnscore_helper import point, atom


class TestPoint(unittest.TestCase):
  """
  Test point.
  """
  def setUp(self):
    """
    Instantiate a local point object.
    """
    self.point_a = point(1,2,3)
    self.point_b = point(-1,-2,-3)

  def testCopyOf(self):
    """
    Verify that copy_of copies x,y,z correctly.
    """
    copy_point = self.point_a.copy_of()
    assert copy_point.x == 1
    assert copy_point.y == 2
    assert copy_point.z == 3

  def testDistTo(self):
    """
    Verify that dist_to implements L2-distance.
    """
    ## || point_a - point_b ||_2 = ||(1,2,3) - (-1,-2,-3)||_2
    ##                           = ||(2,4,6)||_2
    ##                           = sqrt(4 + 16 + 36)
    ##                           = sqrt(56)
    assert self.point_a.dist_to(self.point_b) == np.sqrt(56)

  def testMagnitude(self):
    """
    Verify that magnitude implements L2-Norm.
    """
    ## || (1, 2, 3) ||_2 = || (-1, -2, -3) ||_2
    ##                   = sqrt(1 + 4 + 9)
    ##                   = sqrt(14)
    assert self.point_a.Magnitude() == np.sqrt(14)
    assert self.point_b.Magnitude() == np.sqrt(14)

  def testCreatePDBLine(self):
    """
    Verify that PDB Line is in correct format.
    """
    # TODO(bramsundar): Add a more nontrivial test after looking into
    # PDB standard.
    line = self.point_a.CreatePDBLine(1)
    assert type(line) == str
