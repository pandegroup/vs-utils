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
    TestPDB: Saves dummy PDB to file and verifies that it can be reloaded.
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

  def testConnectedAtomsOfGivenElement(self):
    """
    TestPDB: Verifies that connected atom retrieval works.
    """
    # Verify that no atoms are present when we start.
    assert len(self.pdb.AllAtoms.keys()) == 0
    carbon_atom = atom(element="C")
    oxygen_atom = atom(element="O")
    hydrogen_atom = atom(element="H")

    self.pdb.AddNewAtom(carbon_atom)
    self.pdb.AddNewAtom(oxygen_atom)
    self.pdb.AddNewAtom(hydrogen_atom)

    # We want a carboxyl, so C connects O and H
    carbon_atom.IndicesOfAtomsConnecting = [2,3]
    oxygen_atom.IndicesOfAtomsConnecting = [1]
    hydrogen_atom.IndicesOfAtomsConnecting = [1]

    connected_oxygens = self.pdb.ConnectedAtomsOfGivenElement(1, "O")
    assert len(connected_oxygens) == 1

    connected_hydrogens = self.pdb.ConnectedAtomsOfGivenElement(1, "H")
    assert len(connected_hydrogens) == 1

  def testConnectedHeavyAtoms(self):
    """
    TestPDB: Verifies retrieval of connected heavy atoms.
    """
    # Verify that no atoms are present when we start.
    assert len(self.pdb.AllAtoms.keys()) == 0
    carbon_atom = atom(element="C")
    oxygen_atom = atom(element="O")
    hydrogen_atom = atom(element="H")

    self.pdb.AddNewAtom(carbon_atom)
    self.pdb.AddNewAtom(oxygen_atom)
    self.pdb.AddNewAtom(hydrogen_atom)

    # We want a carboxyl, so C connects O and H
    carbon_atom.IndicesOfAtomsConnecting = [2,3]
    oxygen_atom.IndicesOfAtomsConnecting = [1]
    hydrogen_atom.IndicesOfAtomsConnecting = [1]

    connected_heavy_atoms = self.pdb.ConnectedHeavyAtoms(1)
    assert len(connected_heavy_atoms) == 1
    assert connected_heavy_atoms[0] == 2

  def testCreateNonProteinAtomBondsByDistance(self):
    """
    TestPDB: Verifies creation of bonds.
    """
    carbon_atom = atom(element="C", coordinates=point(0,0,1))
    oxygen_atom = atom(element="O", coordinates=point(0,0,2))

    self.pdb.AddNewNonProteinAtom(carbon_atom)
    self.pdb.AddNewNonProteinAtom(oxygen_atom)

    self.pdb.CreateNonProteinAtomBondsByDistance()
    assert len(carbon_atom.IndicesOfAtomsConnecting) == 1
    assert len(oxygen_atom.IndicesOfAtomsConnecting) == 1

  def testAssignNonProteinCharges(self):
    """
    TestPDB: Verify that non-protein charges are assigned properly.
    """
    # Test metallic ion charge.
    self.pdb = PDB()
    assert len(self.pdb.charges) == 0
    magnesium_atom = atom(element="MG", coordinates=point(0,0,0))
    self.pdb.AddNewNonProteinAtom(magnesium_atom)
    self.pdb.AssignNonProteinCharges()
    assert len(self.pdb.charges) == 1

    # Test ammonium
    self.pdb = PDB()
    assert len(self.pdb.charges) == 0
    # We assign the coordinates to form a tetrahedron; see
    # http://en.wikipedia.org/wiki/Tetrahedron#Formulas_for_a_regular_tetrahedron
    nitrogen_atom = atom(element="N", coordinates=point(0,0,0))
    hydrogen_atom1 = atom(element="H", coordinates=point(1,0,-1./np.sqrt(2)))
    hydrogen_atom2 = atom(element="H", coordinates=point(-1,0,-1./np.sqrt(2)))
    hydrogen_atom3 = atom(element="H", coordinates=point(0,1,1./np.sqrt(2)))
    hydrogen_atom4 = atom(element="H", coordinates=point(0,-1,1./np.sqrt(2)))
    self.pdb.AddNewNonProteinAtom(nitrogen_atom)
    self.pdb.AddNewNonProteinAtom(hydrogen_atom1)
    self.pdb.AddNewNonProteinAtom(hydrogen_atom2)
    self.pdb.AddNewNonProteinAtom(hydrogen_atom3)
    self.pdb.AddNewNonProteinAtom(hydrogen_atom4)
    nitrogen_atom.IndicesOfAtomsConnecting = [2, 3, 4, 5]
    self.pdb.AssignNonProteinCharges()
    assert len(self.pdb.charges) == 1

    # TODO(bramsundar): Figure out how to add more tests here. The issue is
    # that it becomes challenging to specify complicated geometries in the
    # test scripts. Maybe have a collection of data PDBs for tests?

