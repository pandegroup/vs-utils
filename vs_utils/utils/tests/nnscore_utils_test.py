"""
Test NNScore Helper Classes
"""
# pylint mistakenly reports numpy errors:
#     pylint: disable=E1101

import unittest
import numpy as np

from vs_utils.utils.nnscore_utils import Point
from vs_utils.utils.nnscore_utils import Atom
from vs_utils.utils.nnscore_utils import angle_between_three_points
from vs_utils.utils.nnscore_utils import angle_between_points
from vs_utils.utils.nnscore_utils import average_point
from vs_utils.utils.nnscore_utils import cross_product
from vs_utils.utils.nnscore_utils import dot_product
from vs_utils.utils.nnscore_utils import normalized_vector
from vs_utils.utils.nnscore_utils import project_point_onto_plane
from vs_utils.utils.nnscore_utils import vector_scalar_multiply
from vs_utils.utils.nnscore_utils import vector_subtraction


class TestPoint(unittest.TestCase):
  """
  Test point class.
  """
  def setUp(self):
    """
    Instantiate local points for tests.
    """
    self.point_a = Point(coords=np.array([1, 2, 3]))
    self.point_b = Point(coords=np.array([-1, -2, -3]))

  def test_copy_of(self):
    """
    TestPoint: Verify that copy_of copies x,y,z correctly.
    """
    copy_point = self.point_a.copy_of()
    assert np.array_equal(np.array([1, 2, 3]), copy_point.as_array())

  def test_average_point(self):
    """
    TestPoint: Verify that averaging works.
    """
    avg_point = average_point([self.point_a, self.point_b])
    assert np.array_equal(np.array([0, 0, 0]), avg_point.as_array())

  def test_dist_to(self):
    """
    TestPoint: Verify that dist_to implements L2-distance.
    """
    ## || point_a - point_b ||_2 = ||(1,2,3) - (-1,-2,-3)||_2
    ##                           = ||(2,4,6)||_2
    ##                           = sqrt(4 + 16 + 36)
    ##                           = sqrt(56)
    assert self.point_a.dist_to(self.point_b) == np.sqrt(56)

  def test_magnitude(self):
    """
    TestPoint: Verify that magnitude implements L2-Norm.
    """
    ## || (1, 2, 3) ||_2 = || (-1, -2, -3) ||_2
    ##                   = sqrt(1 + 4 + 9)
    ##                   = sqrt(14)
    assert self.point_a.magnitude() == np.sqrt(14)
    assert self.point_b.magnitude() == np.sqrt(14)

  def test_vector_subtraction(self):
    """
    TestPoint: Test that point subtraction works.
    """
    assert np.array_equal(
        vector_subtraction(self.point_a, self.point_b).as_array(),
        np.array([2, 4, 6]))

  def test_cross_product(self):
    """
    TestPoint: Test that cross product works.
    """
    assert np.array_equal(
        cross_product(self.point_a, self.point_b).as_array(),
        np.array([0, 0, 0]))

  def test_vector_scalar(self):
    """
    TestPoint: Test that scalar multiplication works.
    """
    assert np.array_equal(
        vector_scalar_multiply(self.point_a, 2).as_array(),
        np.array([2, 4, 6]))

  def test_dot_prodcut(self):
    """
    TestPoint: Test that dot product works.
    """
    assert dot_product(self.point_a, self.point_b) == -14

  def test_angle_between_three(self):
    """
    TestPoint: Test that the angle between three points is computed.
    """
    assert angle_between_three_points(
        Point(coords=np.array([1, 0, 0])),
        Point(coords=np.array([0, 0, 0])),
        Point(coords=np.array([0, 0, 1]))) == np.pi/2

  def test_angle_between_points(self):
    """
    TestPoint: Test that the angle between two points is computed correctly.
    """
    assert angle_between_points(self.point_a, self.point_b) == np.pi

  def test_normalized_point(self):
    """
    TestPoint: Test that points are normalized.
    """
    assert np.array_equal(
        normalized_vector(Point(coords=np.array([2, 0, 0]))).as_array(),
        np.array([1, 0, 0]))

  def test_project_point(self):
    """
    TestPoint: Test that projection onto plane works.
    """
    # First test with projection onto xy-plane
    value = project_point_onto_plane(
        Point(coords=np.array([1, 2, 3])), [0, 0, 1, 0]) 
    assert np.array_equal(value.as_array(), np.array([1, 2, 0]))

    # Now test projection onto plane z = 4
    value = project_point_onto_plane(
        Point(coords=np.array([1, 2, 3])), [0, 0, 1, 4])
    assert np.array_equal(value.as_array(), np.array([1, 2, 4]))



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
    self.trial_atom.coordinates = Point(coords=np.array([1, 2, 3]))
    self.trial_atom.charge = 0.
    self.trial_atom.element = "C"
    self.trial_atom.residue = "CYS"
    # TODO(bramsundar): Fill in a non-junk value for chain.
    self.trial_atom.chain = "FF"
    self.trial_atom.indices_of_atoms_connecting = [4, 5, 6]

  def test_copy_of(self):
    """
    TestAtom: Verify that copy_of preserves atom information.
    """
    copy_atom = self.trial_atom.copy_of()
    assert copy_atom.atomname == "C"
    assert np.array_equal(copy_atom.coordinates.as_array(), np.array([1, 2, 3]))
    assert copy_atom.charge == 0
    assert copy_atom.element == "C"
    assert copy_atom.residue == "CYS"
    assert copy_atom.chain == "FF"
    assert copy_atom.indices_of_atoms_connecting == [4, 5, 6]

  def test_create_pdb_line(self):
    """
    TestAtom: Verify that PDB Line is in correct format.
    """
    # TODO(bramsundar): Add a more nontrivial test after looking into
    # PDB standard.
    line = self.trial_atom.create_pdb_line(1)
    assert type(line) == str

  def test_number_of_neighors(self):
    """
    TestAtom: Verify that the number of neighbors is computed correctly.
    """
    assert self.empty_atom.number_of_neighbors() == 0
    assert self.trial_atom.number_of_neighbors() == 3

