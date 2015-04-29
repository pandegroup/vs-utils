"""
Test NNScore Helper Classes
"""
import os
import tempfile
import shutil
import unittest
import numpy as np

from vs_utils.features.nnscore_helper import Point, Atom, PDB
from vs_utils.features.tests import __file__ as test_directory


def data_dir():
  """Get location of data directory."""
  return os.path.join(os.path.dirname(test_directory), "data")


class TestPoint(unittest.TestCase):
  """
  Test point class.
  """
  def setUp(self):
    """
    Instantiate local points for tests.
    """
    self.point_a = Point(1,2,3)
    self.point_b = Point(-1,-2,-3)

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
    self.trial_atom.coordinates = Point(1,2,3)
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

    _, self.pdb_filename = tempfile.mkstemp(suffix=".pdb",
        dir=self.temp_dir)

    self.prgr_pdb = PDB()
    prgr_pdb_path = os.path.join(data_dir(), "prgr.pdb")
    self.prgr_pdb.load_PDB_from_file(prgr_pdb_path)

    self.benzene_pdb = PDB()
    benzene_pdb_path = os.path.join(data_dir(), "benzene.pdb")
    self.benzene_pdb.load_PDB_from_file(benzene_pdb_path)

  def tearDown(self):
    """
    Delete temporary directory.
    """
    shutil.rmtree(self.temp_dir)

  def testSaveAndLoad(self):
    """
    TestPDB: Saves dummy PDB to file and verifies that it can be reloaded.
    """
    self.pdb.save_PDB(self.pdb_filename)
    empty_pdb = PDB()
    with open(self.pdb_filename) as pdb_file:
      for line in pdb_file:
        print line
    empty_pdb.load_PDB_from_file(self.pdb_filename)

  def testAddNewAtom(self):
    """
    TestPDB: Verifies that new atoms can be added.
    """
    # Verify that no atoms are present when we start.
    assert len(self.pdb.all_atoms.keys()) == 0
    empty_atom = Atom()
    self.pdb.add_new_atom(empty_atom)
    # Verify that we now have one atom
    assert len(self.pdb.all_atoms.keys()) == 1

  def testConnectedAtomsOfGivenElement(self):
    """
    TestPDB: Verifies that connected atom retrieval works.
    """
    # Verify that no atoms are present when we start.
    assert len(self.pdb.all_atoms.keys()) == 0
    carbon_atom = Atom(element="C")
    oxygen_atom = Atom(element="O")
    hydrogen_atom = Atom(element="H")

    self.pdb.add_new_atom(carbon_atom)
    self.pdb.add_new_atom(oxygen_atom)
    self.pdb.add_new_atom(hydrogen_atom)

    # We want a carboxyl, so C connects O and H
    carbon_atom.indices_of_atoms_connecting = [2,3]
    oxygen_atom.indices_of_atoms_connecting = [1]
    hydrogen_atom.indices_of_atoms_connecting = [1]

    connected_oxygens = self.pdb.connected_atoms_of_given_element(1, "O")
    assert len(connected_oxygens) == 1

    connected_hydrogens = self.pdb.connected_atoms_of_given_element(1, "H")
    assert len(connected_hydrogens) == 1

  def testConnectedHeavyAtoms(self):
    """
    TestPDB: Verifies retrieval of connected heavy atoms.
    """
    # Verify that no atoms are present when we start.
    assert len(self.pdb.all_atoms.keys()) == 0
    carbon_atom = Atom(element="C")
    oxygen_atom = Atom(element="O")
    hydrogen_atom = Atom(element="H")

    self.pdb.add_new_atom(carbon_atom)
    self.pdb.add_new_atom(oxygen_atom)
    self.pdb.add_new_atom(hydrogen_atom)

    # We want a carboxyl, so C connects O and H
    carbon_atom.indices_of_atoms_connecting = [2,3]
    oxygen_atom.indices_of_atoms_connecting = [1]
    hydrogen_atom.indices_of_atoms_connecting = [1]

    connected_heavy_atoms = self.pdb.connected_heavy_atoms(1)
    assert len(connected_heavy_atoms) == 1
    assert connected_heavy_atoms[0] == 2

  def testCreateNonProteinAtomBondsByDistance(self):
    """
    TestPDB: Verifies creation of bonds.
    """
    # First test a toy example
    carbon_atom = Atom(element="C", coordinates=Point(0,0,1))
    oxygen_atom = Atom(element="O", coordinates=Point(0,0,2))

    self.pdb.add_new_non_protein_atom(carbon_atom)
    self.pdb.add_new_non_protein_atom(oxygen_atom)

    self.pdb.create_non_protein_atom_bonds_by_distance()
    assert len(carbon_atom.indices_of_atoms_connecting) == 1
    assert len(oxygen_atom.indices_of_atoms_connecting) == 1

    # Test that all bonds in benzene are created
    assert (len(self.benzene_pdb.all_atoms.keys())
         == len(self.benzene_pdb.non_protein_atoms.keys()))

    for atom_ind in self.benzene_pdb.non_protein_atoms:
      print "Atom %d" % atom_ind
      atom_obj = self.benzene_pdb.non_protein_atoms[atom_ind]
      print "Connected to Atoms: " + str(atom_obj.indices_of_atoms_connecting)

    assert 0 == 1

  def testAssignNonProteinCharges(self):
    """
    TestPDB: Verify that non-protein charges are assigned properly.
    """
    # Test metallic ion charge.
    self.pdb = PDB()
    assert len(self.pdb.charges) == 0
    magnesium_atom = Atom(element="MG", coordinates=Point(0,0,0))
    self.pdb.add_new_non_protein_atom(magnesium_atom)
    self.pdb.assign_non_protein_charges()
    assert len(self.pdb.charges) == 1

    # Test ammonium
    self.pdb = PDB()
    assert len(self.pdb.charges) == 0
    # We assign the coordinates to form a tetrahedron; see
    # http://en.wikipedia.org/wiki/Tetrahedron#Formulas_for_a_regular_tetrahedron
    nitrogen_atom = Atom(element="N", coordinates=Point(0,0,0))
    hydrogen_atom1 = Atom(element="H", coordinates=Point(1,0,-1./np.sqrt(2)))
    hydrogen_atom2 = Atom(element="H", coordinates=Point(-1,0,-1./np.sqrt(2)))
    hydrogen_atom3 = Atom(element="H", coordinates=Point(0,1,1./np.sqrt(2)))
    hydrogen_atom4 = Atom(element="H", coordinates=Point(0,-1,1./np.sqrt(2)))
    self.pdb.add_new_non_protein_atom(nitrogen_atom)
    self.pdb.add_new_non_protein_atom(hydrogen_atom1)
    self.pdb.add_new_non_protein_atom(hydrogen_atom2)
    self.pdb.add_new_non_protein_atom(hydrogen_atom3)
    self.pdb.add_new_non_protein_atom(hydrogen_atom4)
    nitrogen_atom.indices_of_atoms_connecting = [2, 3, 4, 5]
    self.pdb.assign_non_protein_charges()
    assert len(self.pdb.charges) == 1

  def testLigandAssignAromaticRings(self):
    """
    TestPDB: Verify that aromatic rings in ligands are identified.
    """
    # A benzene should have exactly one aromatic ring.
    #print ("len(self.benzene_pdb.aromatic_rings): "
    #  + str(len(self.benzene_pdb.aromatic_rings)))
    #assert len(self.benzene_pdb.aromatic_rings) == 1
