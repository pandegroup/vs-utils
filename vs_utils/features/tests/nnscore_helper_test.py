"""
Test NNScore Helper Classes
"""
import os
import tempfile
import shutil
import unittest
import numpy as np

from vs_utils.features.nnscore_helper import Point
from vs_utils.features.nnscore_helper import Atom
from vs_utils.features.nnscore_helper import PDB
from vs_utils.features.nnscore_helper import average_point
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

  def testAveragePoint(self):
    avg_point = average_point(self.point_a, self.point_b)
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
    assert self.empty_atom.number_of_neighbors() == 0
    assert self.trial_atom.number_of_neighbors() == 3

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

  def testLoadBondsFromPDBList(self):
    """
    TestPDB: Verifies that bonds can be loaded from PDB.
    """
    # Test that we can load CO2
    carbon_atom = Atom(element="C")
    oxygen_atom_1 = Atom(element="O")
    oxygen_atom_2 = Atom(element="O")

    self.pdb.add_new_atom(carbon_atom)
    self.pdb.add_new_atom(oxygen_atom_1)
    self.pdb.add_new_atom(oxygen_atom_2)
    lines = [
      "CONECT    1    2    3                                                 "
      "CONECT    2                                                           "
      "CONECT    3                                                           "
    ]
    self.pdb.load_bonds_from_PDB_list(lines)
    assert len(carbon_atom.indices_of_atoms_connecting) == 2
    assert len(oxygen_atom_1.indices_of_atoms_connecting) == 0
    assert len(oxygen_atom_2.indices_of_atoms_connecting) == 0


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

  def testAssignNonProteinCharges(self):
    """
    TestPDB: Verify that non-protein charges are assigned properly.
    """
    # Test metallic ion charge.
    magnesium_pdb = PDB()
    assert len(magnesium_pdb.charges) == 0
    magnesium_atom = Atom(element="MG", coordinates=Point(0,0,0))
    magnesium_pdb.add_new_non_protein_atom(magnesium_atom)
    magnesium_pdb.assign_non_protein_charges()
    assert len(magnesium_pdb.charges) == 1


  def testIdentifyNitrogenCharges(self):
    """
    TestPDB: Verify that nitrogen groups are charged correctly.
    """
    # Test ammonium sulfate: (NH4)+(NH4)+(SO4)(2-)
    # The labeling should pick up 2 charged nitrogen groups for two
    # ammoniums.
    ammonium_sulfate_pdb = PDB()
    ammonium_sulfate_pdb_path = os.path.join(data_dir(),
        "ammonium_sulfate.pdb")
    ammonium_sulfate_pdb.load_PDB_from_file(
        ammonium_sulfate_pdb_path)
    nitrogen_charges = ammonium_sulfate_pdb.identify_nitrogen_group_charges()
    assert len(nitrogen_charges) == 2
    assert nitrogen_charges[0].positive  # Should be positive
    assert nitrogen_charges[1].positive  # Should be positive

    # Test pyrrolidine (CH2)4NH. The nitrogen here should be sp3
    # hybridized, so is likely to pick up an extra proton to its nitrogen
    # at physiological pH.
    pyrrolidine_pdb = PDB()
    pyrrolidine_pdb_path = os.path.join(data_dir(),
        "pyrrolidine.pdb")
    pyrrolidine_pdb.load_PDB_from_file(pyrrolidine_pdb_path)
    nitrogen_charges = pyrrolidine_pdb.identify_nitrogen_group_charges()
    assert len(nitrogen_charges) == 1
    assert nitrogen_charges[0].positive  # Should be positive

  def testIdentifyCarbonCharges(self):
    """
    TestPDB: Verify that carbon groups are charged correctly.
    """
    # Guanidine is positively charged at physiological pH
    guanidine_pdb = PDB()
    guanidine_pdb_path = os.path.join(data_dir(),
        "guanidine.pdb")
    guanidine_pdb.load_PDB_from_file(
        guanidine_pdb_path)
    carbon_charges = guanidine_pdb.identify_carbon_group_charges()
    assert len(carbon_charges) == 1
    assert carbon_charges[0].positive  # Should be positive

    # sulfaguanidine contains a guanidine group that is likely to be
    # positively protonated at physiological pH
    sulfaguanidine_pdb = PDB()
    sulfaguanidine_pdb_path = os.path.join(data_dir(),
        "sulfaguanidine.pdb")
    sulfaguanidine_pdb.load_PDB_from_file(
        sulfaguanidine_pdb_path)
    carbon_charges = sulfaguanidine_pdb.identify_carbon_group_charges()
    assert len(carbon_charges) == 1
    assert carbon_charges[0].positive  # Should be positive

    # Formic acid is a carboxylic acid, which
    formic_acid_pdb = PDB()
    formic_acid_pdb_path = os.path.join(data_dir(),
        "formic_acid.pdb")
    formic_acid_pdb.load_PDB_from_file(
        formic_acid_pdb_path)
    carbon_charges = formic_acid_pdb.identify_carbon_group_charges()
    assert len(carbon_charges) == 1
    assert not carbon_charges[0].positive  # Should be negatively charged.


  def testLigandAssignAromaticRings(self):
    """
    TestPDB: Verify that aromatic rings in ligands are identified.
    """
    benzene_pdb = PDB()
    benzene_pdb_path = os.path.join(data_dir(), "benzene.pdb")
    benzene_pdb.load_PDB_from_file(benzene_pdb_path)

    # A benzene should have exactly one aromatic ring.
    assert len(benzene_pdb.aromatic_rings) == 1
    # The first 6 atoms in the benzene pdb form the aromatic ring.
    assert (set(benzene_pdb.aromatic_rings[0].indices)
         == set([1,2,3,4,5,6]))
