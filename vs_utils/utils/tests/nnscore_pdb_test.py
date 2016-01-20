"""
Test nnscore pdb implementation. 

Most of the PDB data files for ligands were generated by obabel from the
corresponding pubchem entries. The receptor PDBs (prgr) are taken from the
DUD-E dataset.
"""
# pylint mistakenly reports numpy errors:
#     pylint: disable=E1101

import os
import shutil
import unittest
import numpy as np
import tempfile

from vs_utils.utils.nnscore_pdb import PDB
from vs_utils.utils.nnscore_pdb import remove_redundant_rings 
from vs_utils.utils.nnscore_utils import Point
from vs_utils.utils.nnscore_utils import Atom
from vs_utils.utils.tests import __file__ as test_directory


def data_dir():
  """Get location of data directory."""
  return os.path.join(os.path.dirname(test_directory), "data")


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
    prgr_pdb_path = os.path.join(data_dir(), "prgr_hyd.pdb")
    prgr_pdbqt_path = os.path.join(data_dir(), "prgr_hyd.pdbqt")
    self.prgr_pdb.load_from_files(prgr_pdb_path, prgr_pdbqt_path)

    self._1r5y_protein = PDB()
    _1r5y_protein_pdb = os.path.join(data_dir(), "1r5y_protein_hyd.pdb")
    _1r5y_protein_pdbqt = os.path.join(data_dir(), "1r5y_protein_hyd.pdbqt")
    self._1r5y_protein.load_from_files(_1r5y_protein_pdb, _1r5y_protein_pdbqt)

    self.proteins = [("prgr", self.prgr_pdb), ("1r5y", self._1r5y_protein)]



  def tearDown(self):
    """
    Delete temporary directory.
    """
    shutil.rmtree(self.temp_dir)

  def test_add_new_atom(self):
    """
    TestPDB: Verifies that new atoms can be added.
    """
    # Verify that no atoms are present when we start.
    assert len(self.pdb.all_atoms.keys()) == 0
    empty_atom = Atom()
    self.pdb.add_new_atom(empty_atom)
    # Verify that we now have one atom
    assert len(self.pdb.all_atoms.keys()) == 1

  def test_get_residues(self):
    """
    TestPDB: Tests that all residues in PDB are identified.
    """
    residues = self.prgr_pdb.get_residues()
    # prgr.pdb has 280 unique residues
    assert len(residues.keys()) == 280
    prgr_residues = ["LEU", "ILE", "ASN", "LEU", "LEU", "MET", "SER",
        "ILE", "GLU", "PRO", "ASP", "VAL", "ILE", "TYR", "ALA", "GLY", "HIS",
        "ASP", "THR", "SER", "SER", "SER", "LEU", "LEU", "THR", "SER", "LEU",
        "ASN", "GLN", "LEU", "GLY", "GLU", "ARG", "GLN", "LEU", "LEU", "SER",
        "VAL", "VAL", "LYS", "TRP", "SER", "LYS", "SER", "LEU", "PRO", "GLY",
        "PHE", "ARG", "LEU", "HIS", "ILE", "ASP", "ASP", "GLN", "ILE", "THR",
        "LEU", "ILE", "GLN", "TYR", "SER", "TRP", "MET", "SER", "LEU", "MET",
        "VAL", "PHE", "GLY", "LEU", "GLY", "TRP", "ARG", "SER", "TYR", "LYS",
        "HIS", "VAL", "SER", "GLY", "GLN", "MET", "LEU", "TYR", "PHE", "ALA",
        "PRO", "ASP", "LEU", "ILE", "LEU", "ASN", "GLU", "GLN", "ARG", "MET",
        "LYS", "GLU", "PHE", "TYR", "SER", "LEU", "CYS", "LEU", "THR", "MET",
        "TRP", "GLN", "ILE", "PRO", "GLN", "GLU", "PHE", "VAL", "LYS", "LEU",
        "GLN", "VAL", "SER", "GLN", "GLU", "GLU", "PHE", "LEU", "CYS", "MET",
        "LYS", "VAL", "LEU", "LEU", "LEU", "LEU", "ASN", "THR", "ILE", "PRO",
        "LEU", "GLU", "GLY", "LEU", "PHE", "MET", "ARG", "TYR", "ILE", "GLU",
        "LEU", "ALA", "ILE", "ARG", "ARG", "PHE", "TYR", "GLN", "LEU", "THR",
        "LYS", "LEU", "LEU", "ASP", "ASN", "LEU", "HIS", "ASP", "LEU", "VAL",
        "LYS", "GLN", "LEU", "HIS", "LEU", "TYR", "CYS", "LEU", "ASN", "THR",
        "PHE", "ILE", "GLN", "SER", "ARG", "ALA", "LEU", "SER", "VAL", "GLU",
        "PHE", "PRO", "GLU", "MET", "MET", "SER", "GLU", "VAL", "ILE", "ALA",
        "ALA", "GLN", "LEU", "PRO", "LYS", "ILE", "LEU", "ALA", "GLY", "MET",
        "VAL", "LYS", "PRO", "LEU", "LEU", "PHE", "HIS", "LYS", "ASN", "LEU",
        "ASP", "ASP", "ILE", "THR", "LEU", "ILE", "GLN", "TYR", "SER", "TRP",
        "MET", "THR", "ILE", "PRO", "LEU", "GLU", "GLY", "LEU", "ARG", "VAL",
        "LYS", "GLN", "LEU", "HIS", "LEU", "TYR", "CYS", "LEU", "ASN", "THR",
        "PHE", "ILE", "GLN", "SER", "ARG", "ALA", "LEU", "SER", "VAL", "GLU",
        "PHE", "PRO", "GLU", "MET", "MET", "SER", "GLU", "VAL", "ILE", "ALA",
        "ALA", "GLN", "LEU", "PRO", "LYS", "ILE", "LEU", "ALA", "GLY", "MET",
        "VAL", "LYS", "PRO"]
    # Recall the keys have format RESNAME_RESNUMBER_CHAIN
    resnames = [reskey.split("_")[0].strip() for reskey in residues]
    resnames.sort()
    prgr_residues.sort()
    assert resnames == prgr_residues
    # prgr.pdb has 2749 unique atoms.
    atom_count = 0
    for (_, atom_indices) in residues.iteritems():
      atom_count += len(atom_indices)
    print atom_count
    assert atom_count == 2788

  def test_get_lysine_charges(self):
    """
    TestPDB: Test that lysine charges are identified correctly.
    """
    res_list = self.prgr_pdb.get_residues()
    lysine_charges = self.prgr_pdb.get_lysine_charges(res_list)
    # prgr has 14 lysines.
    print len(lysine_charges)
    assert len(lysine_charges) == 14
    for charge in lysine_charges:
      # Lysine should be posistively charged
      assert charge.positive

  def test_get_arginine_charges(self):
    """
    TestPDB: Test that arginine charges are identified correctly.
    """
    res_list = self.prgr_pdb.get_residues()
    arginine_charges = self.prgr_pdb.get_arginine_charges(res_list)
    # prgr has 10 arginines
    assert len(arginine_charges) == 10
    for charge in arginine_charges:
      # The guanidium in arginine should be positively charged.
      assert charge.positive

  def test_get_histidine_charges(self):
    """
    TestPDB: Test that histidine charges are identified correctly.
    """
    res_list = self.prgr_pdb.get_residues()
    histidine_charges = self.prgr_pdb.get_histidine_charges(res_list)
    # prgr has 7 arginines
    assert len(histidine_charges) == 7
    for charge in histidine_charges:
      # The nitrogens pick up positive charges
      assert charge.positive

  def test_get_glutamic_acid_charges(self):
    """
    TestPDB: Test that glutamic acid charges are identified correctly.
    """
    res_list = self.prgr_pdb.get_residues()
    glutamic_acid_charges = self.prgr_pdb.get_glutamic_acid_charges(res_list)
    assert len(glutamic_acid_charges) == 16
    for charge in glutamic_acid_charges:
      # The carboxyls get deprotonated.
      assert not charge.positive

  def test_get_aspartic_acid_charges(self):
    """
    TestPDB: Test that aspartic acid charges are identified correctly.
    """
    res_list = self.prgr_pdb.get_residues()
    aspartic_acid_charges = self.prgr_pdb.get_aspartic_acid_charges(res_list)
    assert len(aspartic_acid_charges) == 9
    for charge in aspartic_acid_charges:
      # The carboxyls get deprotonated
      assert not charge.positive

  def test_assign_ligand_aromatics(self):
    """
    TestPDB: Test that non-protein aromatic rings are assigned correctly.
    """
    ### 3ao4 comes from PDBBind-CN and contains some cruft in the PDB file:
    ### atoms without residues labelled. This triggered some problems with
    ### non-protein aromatics complaining.
    # TODO(rbharath): Add a stub here.
    _3ao4_protein = PDB()
    _3ao4_protein_pdb = os.path.join(data_dir(), "3ao4_protein_hyd.pdb")
    _3ao4_protein_pdbqt = os.path.join(data_dir(), "3ao4_protein_hyd.pdbqt")
    _3ao4_protein.load_from_files(_3ao4_protein_pdb, _3ao4_protein_pdbqt)

  def test_remove_redundant_rings(self):
    """
    TestPDB: Test that redundant rings are removed.
    """
    # Recall that each ring is represented as a list of atom indices.
    # Test that rings of length 0 are removed
    assert remove_redundant_rings([[]]) == []
    # Set that supersets are removed
    assert (remove_redundant_rings([[1, 2, 3], [1, 3, 4, 5], [1, 2, 3, 4, 5]])
        == [[1, 2, 3], [1, 3, 4, 5]])
    # Ensure that duplicate rings are handled correctly (that is, only one
    # copy of a duplicate ring should remain)
    assert remove_redundant_rings([[1, 2, 3], [1, 3, 2]]) == [[1, 2, 3]]

  def test_assign_protein_aromatics(self):
    """
    TestPDB: Test that aromatic rings are assigned correctly.
    """
    for name, protein in self.proteins:
      # The proteins should have aromatic rings assigned already by
      # load_from_files()
      print "Processing aromatics for %s" % name
      for aromatic in protein.aromatic_rings:
        assert aromatic is not None

  def test_get_phenylalanine_aromatics(self):
    """
    TestPDB: Test that phenylalanine aromatic rings are retrieved.
    """
    res_list = self.prgr_pdb.get_residues()
    phenylalanine_aromatics = (
        self.prgr_pdb.get_phenylalanine_aromatics(res_list))

    # prgr has 13 phenylalanines, each of which has 1 aromatic ring.
    assert len(phenylalanine_aromatics) == 13
    for aromatic in phenylalanine_aromatics:
      # The aromatic rings in phenylalanine have 6 elements each
      assert len(aromatic.indices) == 6

  def test_get_tyrosine_aromatics(self):
    """
    TestPDB: Test that tyrosine aromatic rings are retrieved.
    """
    # prgr has 10 tyrosines, each of which has 1 aromatic ring.
    res_list = self.prgr_pdb.get_residues()
    tyrosine_aromatics = self.prgr_pdb.get_tyrosine_aromatics(res_list)
    assert len(tyrosine_aromatics) == 10
    for aromatic in tyrosine_aromatics:
      # The aromatic rings in tyrosine have 6 elements each
      assert len(aromatic.indices) == 6

  def test_get_histidine_aromatics(self):
    """
    TestPDB: Test that histidine aromatic rings are retrieved.
    """
    res_list = self.prgr_pdb.get_residues()
    histidine_aromatics = self.prgr_pdb.get_histidine_aromatics(res_list)
    # prgr has 7 histidines, each of which has 1 aromatic ring.
    assert len(histidine_aromatics) == 7
    for aromatic in histidine_aromatics:
      # The aromatic rings in histidine have 6 elements each
      print len(aromatic.indices)
      assert len(aromatic.indices) == 5

  def test_get_tryptophan_aromatics(self):
    """
    TestPDB: Test that tryptophan aromatic rings are retrieved.
    """
    res_list = self.prgr_pdb.get_residues()
    tryptophan_aromatics = self.prgr_pdb.get_tryptophan_aromatics(res_list)
    # prgr has 5 tryptophans, each of which has 2 aromatic ring.
    print len(tryptophan_aromatics)
    assert len(tryptophan_aromatics) == 10 
    num_five_rings, num_six_rings = 0, 0
    for aromatic in tryptophan_aromatics:
      # One aromatic ring in tryptophan hahas 6 elements each,
      # while the other has 5 elements.
      if len(aromatic.indices) == 6:
        num_six_rings += 1
      elif len(aromatic.indices) == 5:
        num_five_rings += 1
    assert num_six_rings == 5
    assert num_five_rings == 5

  def test_connected_atoms(self):
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
    carbon_atom.indices_of_atoms_connecting = [2, 3]
    oxygen_atom.indices_of_atoms_connecting = [1]
    hydrogen_atom.indices_of_atoms_connecting = [1]

    connected_oxygens = self.pdb.connected_atoms(1, "O")
    assert len(connected_oxygens) == 1

    connected_hydrogens = self.pdb.connected_atoms(1, "H")
    assert len(connected_hydrogens) == 1

  def test_load_bonds_from_pdb(self):
    """
    TestPDB: Verifies that bonds can be loaded from PDB.
    """
    pdb = PDB()
    # Test that we can load CO2
    carbon_atom = Atom(element="C")
    oxygen_atom_1 = Atom(element="O")
    oxygen_atom_2 = Atom(element="O")

    pdb.add_new_atom(carbon_atom)
    pdb.add_new_atom(oxygen_atom_1)
    pdb.add_new_atom(oxygen_atom_2)
    lines = [
      "CONECT    1    2    3                                                 "
      "CONECT    2                                                           "
      "CONECT    3                                                           "
    ]
    with tempfile.NamedTemporaryFile() as temp:
      temp.write("\n".join(lines))
      temp.flush()
      pdb.load_bonds_from_pdb(temp.name)
    assert len(carbon_atom.indices_of_atoms_connecting) == 2
    assert len(oxygen_atom_1.indices_of_atoms_connecting) == 0
    assert len(oxygen_atom_2.indices_of_atoms_connecting) == 0


  def test_connected_heavy_atoms(self):
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
    carbon_atom.indices_of_atoms_connecting = [2, 3]
    oxygen_atom.indices_of_atoms_connecting = [1]
    hydrogen_atom.indices_of_atoms_connecting = [1]

    connected_heavy_atoms = self.pdb.connected_heavy_atoms(1)
    assert len(connected_heavy_atoms) == 1
    assert connected_heavy_atoms[0] == 2

  def test_assign_non_protein_charges(self):
    """
    TestPDB: Verify that charges are properly added to ligands.
    """
    # Test ammonium sulfate: (NH4)+(NH4)+(SO4)(2-)
    # There should be 3 charged groups, two positive, one negative
    ammonium_sulfate_pdb = PDB()
    ammonium_sulfate_pdb_path = os.path.join(data_dir(),
        "ammonium_sulfate_hyd.pdb")
    ammonium_sulfate_pdbqt_path = os.path.join(data_dir(),
        "ammonium_sulfate_hyd.pdbqt")
    # Notice that load automatically identifies non-protein charges.
    ammonium_sulfate_pdb.load_from_files(
        ammonium_sulfate_pdb_path, ammonium_sulfate_pdbqt_path)
    assert len(ammonium_sulfate_pdb.charges) == 3
    num_pos, num_neg = 0, 0
    for charge in ammonium_sulfate_pdb.charges:
      if charge.positive:
        num_pos += 1
      else:
        num_neg += 1
    assert num_pos == 2
    assert num_neg == 1

  def test_metallic_charges(self):
    """
    TestPDB: Verify that non-protein charges are assigned properly.
    """
    # Test metallic ion charge.
    magnesium_pdb = PDB()
    magnesium_atom = Atom(element="MG",
        coordinates=Point(coords=np.array([0,0,0])))
    magnesium_pdb.add_new_non_protein_atom(magnesium_atom)
    metallic_charges = magnesium_pdb.identify_metallic_charges()
    assert len(metallic_charges) == 1

  def test_nitrogen_charges(self):
    """
    TestPDB: Verify that nitrogen groups are charged correctly.
    """
    # Test ammonium sulfate: (NH4)+(NH4)+(SO4)(2-)
    # The labeling should pick up 2 charged nitrogen groups for two
    # ammoniums.
    ammonium_sulfate_pdb = PDB()
    ammonium_sulfate_pdb_path = os.path.join(data_dir(),
        "ammonium_sulfate_hyd.pdb")
    ammonium_sulfate_pdbqt_path = os.path.join(data_dir(),
        "ammonium_sulfate_hyd.pdbqt")
    ammonium_sulfate_pdb.load_from_files(
        ammonium_sulfate_pdb_path, ammonium_sulfate_pdbqt_path)
    nitrogen_charges = ammonium_sulfate_pdb.identify_nitrogen_charges()
    assert len(nitrogen_charges) == 2
    assert nitrogen_charges[0].positive  # Should be positive
    assert nitrogen_charges[1].positive  # Should be positive

    # Test pyrrolidine (CH2)4NH. The nitrogen here should be sp3
    # hybridized, so is likely to pick up an extra proton to its nitrogen
    # at physiological pH.
    pyrrolidine_pdb = PDB()
    pyrrolidine_pdb_path = os.path.join(data_dir(),
        "pyrrolidine_hyd.pdb")
    pyrrolidine_pdbqt_path = os.path.join(data_dir(),
        "pyrrolidine_hyd.pdbqt")
    pyrrolidine_pdb.load_from_files(pyrrolidine_pdb_path,
        pyrrolidine_pdbqt_path)
    nitrogen_charges = pyrrolidine_pdb.identify_nitrogen_charges()
    assert len(nitrogen_charges) == 1
    assert nitrogen_charges[0].positive  # Should be positive

  def test_carbon_charges(self):
    """
    TestPDB: Verify that carbon groups are charged correctly.
    """
    # Guanidine is positively charged at physiological pH
    guanidine_pdb = PDB()
    guanidine_pdb_path = os.path.join(data_dir(),
        "guanidine_hyd.pdb")
    guanidine_pdbqt_path = os.path.join(data_dir(),
        "guanidine_hyd.pdbqt")
    guanidine_pdb.load_from_files(
        guanidine_pdb_path, guanidine_pdbqt_path)
    carbon_charges = guanidine_pdb.identify_carbon_charges()
    assert len(carbon_charges) == 1
    assert carbon_charges[0].positive  # Should be positive

    # sulfaguanidine contains a guanidine group that is likely to be
    # positively protonated at physiological pH
    sulfaguanidine_pdb = PDB()
    sulfaguanidine_pdb_path = os.path.join(data_dir(),
        "sulfaguanidine_hyd.pdb")
    sulfaguanidine_pdbqt_path = os.path.join(data_dir(),
        "sulfaguanidine_hyd.pdbqt")
    sulfaguanidine_pdb.load_from_files(
        sulfaguanidine_pdb_path, sulfaguanidine_pdbqt_path)
    carbon_charges = sulfaguanidine_pdb.identify_carbon_charges()
    assert len(carbon_charges) == 1
    assert carbon_charges[0].positive  # Should be positive

    # Formic acid is a carboxylic acid, which should be negatively charged.
    formic_acid_pdb = PDB()
    formic_acid_pdb_path = os.path.join(data_dir(),
        "formic_acid_hyd.pdb")
    formic_acid_pdbqt_path = os.path.join(data_dir(),
        "formic_acid_hyd.pdbqt")
    formic_acid_pdb.load_from_files(
        formic_acid_pdb_path, formic_acid_pdbqt_path)
    carbon_charges = formic_acid_pdb.identify_carbon_charges()
    assert len(carbon_charges) == 1
    assert not carbon_charges[0].positive  # Should be negatively charged.

  def test_phosphorus_charges(self):
    """
    TestPDB: Verify that Phosphorus groups are charged correctly.
    """
    # CID82671 contains a phosphate between two aromatic groups.
    phosphate_pdb = PDB()
    phosphate_pdb_path = os.path.join(data_dir(),
      "82671_hyd.pdb")
    phosphate_pdbqt_path = os.path.join(data_dir(),
      "82671_hyd.pdb")
    phosphate_pdb.load_from_files(
        phosphate_pdb_path, phosphate_pdbqt_path)
    phosphorus_charges = phosphate_pdb.identify_phosphorus_charges()
    assert len(phosphorus_charges) == 1
    assert not phosphorus_charges[0].positive  # Should be negatively charged.


  def test_sulfur_charges(self):
    """
    TestPDB: Verify that sulfur groups are charged correctly.
    """
    triflic_acid_pdb = PDB()
    triflic_acid_pdb_path = os.path.join(data_dir(),
      "triflic_acid_hyd.pdb")
    triflic_acid_pdbqt_path = os.path.join(data_dir(),
      "triflic_acid_hyd.pdbqt")
    triflic_acid_pdb.load_from_files(
      triflic_acid_pdb_path,
      triflic_acid_pdbqt_path)
    sulfur_charges = (
        triflic_acid_pdb.identify_sulfur_charges())
    assert len(sulfur_charges) == 1
    assert not sulfur_charges[0].positive  # Should be negatively charged.


  def test_ligand_assign_aromatics(self):
    """
    TestPDB: Verify that aromatic rings in ligands are identified.
    """
    benzene_pdb = PDB()
    benzene_pdb_path = os.path.join(data_dir(), "benzene_hyd.pdb")
    benzene_pdbqt_path = os.path.join(data_dir(), "benzene_hyd.pdbqt")
    benzene_pdb.load_from_files(benzene_pdb_path, benzene_pdbqt_path)

    # A benzene should have exactly one aromatic ring.
    print benzene_pdb.aromatic_rings
    assert len(benzene_pdb.aromatic_rings) == 1
    # The first 6 atoms in the benzene pdb form the aromatic ring.
    assert (set(benzene_pdb.aromatic_rings[0].indices)
         == set([1,2,3,4,5,6]))

  def test_assign_secondary_structure(self):
    """
    TestPDB: Verify that secondary structure is assigned meaningfully.
    """
    # TODO(rbharath): This test is just a stub. Add a more realistic test
    # that checks that nontrivial secondary structure is computed correctly
    # here.
    self.prgr_pdb.assign_secondary_structure()
    

  def test_get_structure_dict(self):
    """
    TestPDB: Verify that dict with rudimentary structure labels is generated.

    TODO(rbharath): This is just a stub. Add some nontrivial tests here.
    """
    structures = self.prgr_pdb.get_structure_dict()
    print structures
    print len(structures)