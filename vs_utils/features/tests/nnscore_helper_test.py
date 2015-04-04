"""
Test NNScore Helper Classes
"""
import tempfile
import shutil
import unittest
import numpy as np

from vs_utils.features.nnscore_helper import point, atom, PDB


class TestPoint(unittest.TestCase):
  """
  Test point class.
  """
  def setUp(self):
    """
    Instantiate local points for tests.
    """
    self.point_a = point(1,2,3)
    self.point_b = point(-1,-2,-3)

  def testCopyOf(self):
    """
    TestPoint: Verify that copy_of copies x,y,z correctly.
    """
    copy_point = self.point_a.copy_of()
    assert copy_point.x == 1
    assert copy_point.y == 2
    assert copy_point.z == 3

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
    assert self.point_a.Magnitude() == np.sqrt(14)
    assert self.point_b.Magnitude() == np.sqrt(14)

  def testCreatePDBLine(self):
    """
    TestPoint: Verify that PDB Line is in correct format.
    """
    # TODO(bramsundar): Add a more nontrivial test after looking into
    # PDB standard.
    line = self.point_a.CreatePDBLine(1)
    assert type(line) == str

class TestAtom(unittest.TestCase):
  """
  Test atom class.
  """

  def setUp(self):
    """
    Instantiates a pair of atom objects for tests.
    """
    self.empty_atom = atom()
    self.trial_atom = atom()
    self.trial_atom.atomname = "C"
    self.trial_atom.coordinates = point(1,2,3)
    self.trial_atom.charge = 0.
    self.trial_atom.element = "C"
    self.trial_atom.residue = "CYS"
    # TODO(bramsundar): Fill in a non-junk value for chain.
    self.trial_atom.chain = "FF"
    self.trial_atom.IndicesOfAtomsConnecting = [4, 5, 6]

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
    assert copy_atom.IndicesOfAtomsConnecting == [4, 5, 6]

  def testCreatePDBLine(self):
    """
    TestAtom: Verify that PDB Line is in correct format.
    """
    # TODO(bramsundar): Add a more nontrivial test after looking into
    # PDB standard.
    line = self.trial_atom.CreatePDBLine(1)
    assert type(line) == str

  def testNumberOfNeighors(self):
    """
    TestAtom: Verify that the number of neighbors is computed correctly.
    """
    assert self.empty_atom.NumberOfNeighbors() == 0
    assert self.trial_atom.NumberOfNeighbors() == 3

  # TODO(bramsundar): Add more tests here.

class TestPDB(unittest.TestCase):
  """"
  Test PDB class.
  """

  def setUp(self):
    """
    Instantiate a dummy PDB file.
    """
    self.temp_dir = tempfile.mkdtemp()
    self.pdb = PDB()

    _, self.pdb_filename = tempfile.mkstemp(suffix='.pdb', 
        dir=self.temp_dir)

  def tearDown(self):
    """
    Delete temporary directory.
    """
    shutil.rmtree(self.temp_dir)

  def testSaveAndLoad(self):
    """
    TestPDB: Saves dummpy PDB to file and verifies that it can be reloaded.
    """
    self.pdb.SavePDB(self.pdb_filename) 
    empty_pdb = PDB()
    with open(self.pdb_filename) as pdb_file:
      for line in pdb_file:
        print line
    empty_pdb.LoadPDB_from_file(self.pdb_filename)

  def testAddNewAtom(self):
    """
    TestPDB: Verifies that new atoms can be added.
    """
    # Verify that no atoms are present when we start.
    assert len(self.pdb.AllAtoms.keys()) == 0
    empty_atom = atom()
    self.pdb.AddNewAtom(empty_atom)
    # Verify that we now have one atom
    assert len(self.pdb.AllAtoms.keys()) == 1

