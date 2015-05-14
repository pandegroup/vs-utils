"""
Test NNScore Helper Classes
"""
import os
import tempfile
import shutil
import unittest
import numpy as np

from vs_utils.utils.nnscore_pdb import PDB
from vs_utils.utils.nnscore_utils import Point
from vs_utils.utils.nnscore_utils import Atom
from vs_utils.utils.nnscore_utils import average_point
from vs_utils.utils.tests import __file__ as test_directory


class TestPoint(unittest.TestCase):
  """
  Test point class.
  """
  def setUp(self):
    """
    Instantiate local points for tests.
    """
    self.point_a = Point(coords=np.array([1,2,3]))
    self.point_b = Point(coords=np.array([-1,-2,-3]))

  def testCopyOf(self):
    """
    TestPoint: Verify that copy_of copies x,y,z correctly.
    """
    copy_point = self.point_a.copy_of()
    assert copy_point.x == 1
    assert copy_point.y == 2
    assert copy_point.z == 3

  def testAveragePoint(self):
    avg_point = average_point([self.point_a, self.point_b])
    assert avg_point.x == 0
    assert avg_point.y == 0
    assert avg_point.z == 0

  def testDistTo(self):
    """
    TestPoint: Verify that dist_to implements L2-distance.
    """
    ## || point_a - point_b ||_2 = ||(1,2,3) - (-1,-2,-3)||_2
    ##                           = ||(2,4,6)||_2
    ##                           = sqrt(4 + 16 + 36)
    ##                           = sqrt(56)
    assert self.point_a.dist_to(self.point_b) == np.sqrt(56)

  def testMagnitude(self):
    """
    TestPoint: Verify that magnitude implements L2-Norm.
    """
    ## || (1, 2, 3) ||_2 = || (-1, -2, -3) ||_2
    ##                   = sqrt(1 + 4 + 9)
    ##                   = sqrt(14)
    assert self.point_a.magnitude() == np.sqrt(14)
    assert self.point_b.magnitude() == np.sqrt(14)

class TestAtom(unittest.TestCase):
  """
  Test atom class.
  """

  def setUp(self):
    """
    Instantiates a pair of atom objects for tests.
    """
    self.empty_atom = Atom()
    self.trial_atom = Atom()
    self.trial_atom.atomname = "C"
    self.trial_atom.coordinates = Point(coords=np.array([1,2,3]))
    self.trial_atom.charge = 0.
    self.trial_atom.element = "C"
    self.trial_atom.residue = "CYS"
    # TODO(bramsundar): Fill in a non-junk value for chain.
    self.trial_atom.chain = "FF"
    self.trial_atom.indices_of_atoms_connecting = [4, 5, 6]

  def testCopyOf(self):
    """
    TestAtom: Verify that copy_of preserves atom information.
    """
    copy_atom = self.trial_atom.copy_of()
    assert copy_atom.atomname == "C"
    assert copy_atom.coordinates.x == 1
    assert copy_atom.coordinates.y == 2
    assert copy_atom.coordinates.z == 3
    assert copy_atom.charge == 0
    assert copy_atom.element == "C"
    assert copy_atom.residue == "CYS"
    assert copy_atom.chain == "FF"
    assert copy_atom.indices_of_atoms_connecting == [4, 5, 6]

  def testCreatePDBLine(self):
    """
    TestAtom: Verify that PDB Line is in correct format.
    """
    # TODO(bramsundar): Add a more nontrivial test after looking into
    # PDB standard.
    line = self.trial_atom.create_PDB_line(1)
    assert type(line) == str

  def testNumberOfNeighors(self):
    """
    TestAtom: Verify that the number of neighbors is computed correctly.
    """
    assert self.empty_atom.number_of_neighbors() == 0
    assert self.trial_atom.number_of_neighbors() == 3

