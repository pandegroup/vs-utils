"""
Test NNScore Binana featurizer.

TODO(rbharath): There still isn't an example structure that exhibits
salt-bridge interactions. There might be a bug in the pi-T interaction
finger, and the H-bonds are known to miss some potential bonds with an
overly-conservative bond-angle cutoff. 
"""
# pylint mistakenly reports numpy errors:
#     pylint: disable=E1101
import os
import numpy as np
import unittest

from vs_utils.features.nnscore import Binana
from vs_utils.features.nnscore import compute_hydrophobic_contacts 
from vs_utils.features.nnscore import compute_electrostatic_energy
from vs_utils.features.nnscore import compute_ligand_atom_counts
from vs_utils.features.nnscore import compute_active_site_flexibility
from vs_utils.features.nnscore import compute_pi_t
from vs_utils.features.nnscore import compute_pi_cation
from vs_utils.features.nnscore import compute_pi_pi_stacking
from vs_utils.features.nnscore import compute_hydrogen_bonds
from vs_utils.features.nnscore import compute_contacts
from vs_utils.features.nnscore import compute_salt_bridges
from vs_utils.utils.nnscore_pdb import PDB
from vs_utils.utils.tests import __file__ as test_directory

def data_dir():
  """Get location of data directory."""
  return os.path.join(os.path.dirname(test_directory), "data")

class TestBinana(unittest.TestCase):
  """
  Test Binana Binding Pose Featurizer.
  """
  def setUp(self):
    """
    Instantiate local copy of Binana object.
    """
    self.binana = Binana()

    ### PRGR is from the DUD-E collection
    prgr_receptor = PDB()
    prgr_pdb = os.path.join(data_dir(), "prgr_hyd.pdb")
    prgr_pdbqt = os.path.join(data_dir(), "prgr_hyd.pdbqt")
    prgr_receptor.load_from_files(prgr_pdb, prgr_pdbqt)
    # This compound is CHEMBL1164248
    prgr_active = PDB()
    prgr_active_pdb = os.path.join(data_dir(), "prgr_active0_hyd.pdb")
    prgr_active_pdbqt = os.path.join(data_dir(), "prgr_active0_hyd.pdbqt")
    prgr_active.load_from_files(prgr_active_pdb, prgr_active_pdbqt)

    ### c-Abl is taken from the Autodock Vina examples
    cabl_receptor = PDB()
    cabl_receptor_pdb = os.path.join(data_dir(), "c-Abl_hyd.pdb")
    cabl_receptor_pdbqt = os.path.join(data_dir(), "c-Abl_hyd.pdbqt")
    cabl_receptor.load_from_files(cabl_receptor_pdb,
        cabl_receptor_pdbqt)
    # This compound is imatinib
    cabl_active = PDB()
    cabl_active_pdb = os.path.join(data_dir(), "imatinib_hyd.pdb")
    cabl_active_pdbqt = os.path.join(data_dir(), "imatinib_hyd.pdbqt")
    cabl_active.load_from_files(cabl_active_pdb,
        cabl_active_pdbqt)

    ### 1zea comes from PDBBind-CN
    # Python complains about variables starting with numbers, so put an
    # underscore in front of everything.
    _1zea_protein = PDB()
    _1zea_protein_pdb = os.path.join(data_dir(), "1zea_protein_hyd.pdb")
    _1zea_protein_pdbqt = os.path.join(data_dir(), "1zea_protein_hyd.pdbqt")
    _1zea_protein.load_from_files(_1zea_protein_pdb, _1zea_protein_pdbqt)
    # The ligand is also specified by pdbbind
    _1zea_ligand = PDB()
    _1zea_ligand_pdb = os.path.join(data_dir(), "1zea_ligand_hyd.pdb")
    _1zea_ligand_pdbqt = os.path.join(data_dir(), "1zea_ligand_hyd.pdbqt")
    _1zea_ligand.load_from_files(_1zea_ligand_pdb, _1zea_ligand_pdbqt)

    ### 1r5y comes from PDBBind-CN
    _1r5y_protein = PDB()
    _1r5y_protein_pdb = os.path.join(data_dir(), "1r5y_protein_hyd.pdb")
    _1r5y_protein_pdbqt = os.path.join(data_dir(), "1r5y_protein_hyd.pdbqt")
    _1r5y_protein.load_from_files(_1r5y_protein_pdb, _1r5y_protein_pdbqt)
    # The ligand is also specified by pdbbind
    _1r5y_ligand = PDB()
    _1r5y_ligand_pdb = os.path.join(data_dir(), "1r5y_ligand_hyd.pdb")
    _1r5y_ligand_pdbqt = os.path.join(data_dir(), "1r5y_ligand_hyd.pdbqt")
    _1r5y_ligand.load_from_files(_1r5y_ligand_pdb, _1r5y_ligand_pdbqt)

    ### 3ao4 comes from PDBBind-CN
    _3ao4_protein = PDB()
    _3ao4_protein_pdb = os.path.join(data_dir(), "3ao4_protein_hyd.pdb")
    _3ao4_protein_pdbqt = os.path.join(data_dir(), "3ao4_protein_hyd.pdbqt")
    _3ao4_protein.load_from_files(_3ao4_protein_pdb, _3ao4_protein_pdbqt)
    # The ligand is also specified by pdbbind
    _3ao4_ligand = PDB()
    _3ao4_ligand_pdb = os.path.join(data_dir(), "3ao4_ligand_hyd.pdb")
    _3ao4_ligand_pdbqt = os.path.join(data_dir(), "3ao4_ligand_hyd.pdbqt")
    _3ao4_ligand.load_from_files(_3ao4_ligand_pdb, _3ao4_ligand_pdbqt)

    ### 2jdm comes from PDBBind-CN
    _2jdm_protein = PDB()
    _2jdm_protein_pdb = os.path.join(data_dir(), "2jdm_protein_hyd.pdb")
    _2jdm_protein_pdbqt = os.path.join(data_dir(), "2jdm_protein_hyd.pdbqt")
    _2jdm_protein.load_from_files(_2jdm_protein_pdb, _2jdm_protein_pdbqt)
    # The ligand is also specified by pdbbind
    _2jdm_ligand = PDB()
    _2jdm_ligand_pdb = os.path.join(data_dir(), "2jdm_ligand_hyd.pdb")
    _2jdm_ligand_pdbqt = os.path.join(data_dir(), "2jdm_ligand_hyd.pdbqt")
    _2jdm_ligand.load_from_files(_2jdm_ligand_pdb, _2jdm_ligand_pdbqt)


    self.test_cases = [("prgr", prgr_receptor, prgr_active),
                       ("cabl", cabl_receptor, cabl_active),
                       ("1zea", _1zea_protein, _1zea_ligand),
                       ("1r5y", _1r5y_protein, _1r5y_ligand),
                       ("3ao4", _3ao4_protein, _3ao4_ligand),
                       ("2jdm", _2jdm_protein, _2jdm_ligand)]

  def test_compute_hydrophobic(self):
    """
    TestBinana: Test that hydrophobic contacts are established.
    """
    hydrophobics_dict = {}
    for name, protein, ligand in self.test_cases:
      hydrophobics_dict[name] = compute_hydrophobic_contacts(
          ligand, protein)
    for name, hydrophobics in hydrophobics_dict.iteritems():
      print "Processing hydrohobics for %s" % name
      assert len(hydrophobics) == 6
      assert "BACKBONE_ALPHA" in hydrophobics
      assert "BACKBONE_BETA" in hydrophobics
      assert "BACKBONE_OTHER" in hydrophobics
      assert "SIDECHAIN_ALPHA" in hydrophobics
      assert "SIDECHAIN_BETA" in hydrophobics
      assert "SIDECHAIN_OTHER" in hydrophobics

  def test_compute_electrostatics(self):
    """
    TestBinana: Test that electrostatic energies are computed.
    """
    electrostatics_dict = {}
    for name, protein, ligand in self.test_cases:
      electrostatics_dict[name] = compute_electrostatic_energy(
          ligand, protein)
    for name, electrostatics in electrostatics_dict.iteritems():
      print "Processing electrostatics for %s" % name
      # The keys of these dicts are pairs of atomtypes, but the keys are
      # sorted so that ("C", "O") is always written as "C_O". Thus, for N
      # atom types, there are N*(N+1)/2 unique pairs.
      num_atoms = len(Binana.atom_types)
      assert len(electrostatics) == num_atoms*(num_atoms+1)/2
      assert np.count_nonzero(np.array(electrostatics.values())) > 0

  def test_compute_flexibility(self):
    """
    TestBinana: Gather statistics about active site protein atoms.
    """
    active_site_dict = {}
    for name, protein, ligand in self.test_cases:
      active_site_dict[name] = compute_active_site_flexibility(
          ligand, protein)
    for name, active_site_flexibility in active_site_dict.iteritems():
      print "Processing active site flexibility for %s" % name
      assert len(active_site_flexibility.keys()) == 6
      assert "BACKBONE_ALPHA" in active_site_flexibility
      assert "BACKBONE_BETA" in active_site_flexibility
      assert "BACKBONE_OTHER" in active_site_flexibility
      assert "SIDECHAIN_ALPHA" in active_site_flexibility
      assert "SIDECHAIN_BETA" in active_site_flexibility
      assert "SIDECHAIN_OTHER" in active_site_flexibility

  def test_compute_hydrogen_bonds(self):
    """
    TestBinana: Compute the number of hydrogen bonds.

    TODO(rbharath): The hydrogen-bond angle cutoff seems like it's
    incorrect to me. The hydrogens are placed by openbabel and aren't
    optimized, so I'm pretty sure that this code will miss many hydrogen
    bonds.
    Here are some options:
    -) Find a method to optimize the hydrogen placement.
    -) Place a more permissive angle cutoff for hydrogens.
    -) Allow for "buckets": angles 0-20, 20-40, 40-60, etc. and count the
    number of hydrogen bonds in each bucket.
    """
    hbonds_dict = {}
    for name, protein, ligand in self.test_cases:
      hbonds_dict[name] = compute_hydrogen_bonds(
          ligand, protein)
    for name, hbonds in hbonds_dict.iteritems():
      print "Processing hydrogen bonds for %s" % name
      assert len(hbonds) == 12
      assert "HDONOR-LIGAND_BACKBONE_ALPHA" in hbonds
      assert "HDONOR-LIGAND_BACKBONE_BETA" in hbonds
      assert "HDONOR-LIGAND_BACKBONE_OTHER" in hbonds
      assert "HDONOR-LIGAND_SIDECHAIN_ALPHA" in hbonds
      assert "HDONOR-LIGAND_SIDECHAIN_BETA" in hbonds
      assert "HDONOR-LIGAND_SIDECHAIN_OTHER" in hbonds
      assert "HDONOR-RECEPTOR_BACKBONE_ALPHA" in hbonds
      assert "HDONOR-RECEPTOR_BACKBONE_BETA" in hbonds
      assert "HDONOR-RECEPTOR_BACKBONE_OTHER" in hbonds
      assert "HDONOR-RECEPTOR_SIDECHAIN_ALPHA" in hbonds
      assert "HDONOR-RECEPTOR_SIDECHAIN_BETA" in hbonds
      assert "HDONOR-RECEPTOR_SIDECHAIN_OTHER" in hbonds

  def test_compute_ligand_atom_counts(self):
    """
    TestBinana: Compute the Number of Ligand Atom Counts.
    """
    counts_dict = {}
    for name, _, ligand in self.test_cases:
      counts_dict[name] = compute_ligand_atom_counts(
          ligand)
    for name, counts in counts_dict.iteritems():
      print "Processing ligand atom counts for %s" % name
      assert len(counts) == len(Binana.atom_types)

  def test_compute_contacts(self):
    """
    TestBinana: Compute contacts between Ligand and receptor.
    """
    contacts_dict = {}
    for name, protein, ligand in self.test_cases:
      contacts_dict[name] = compute_contacts(
          ligand, protein)
    for name, (close_contacts, contacts) in contacts_dict.iteritems():
      print "Processing contacts for %s" % name
      print "close_contacts"
      for key, val in close_contacts.iteritems():
        if val != 0:
          print (key, val)
      print "len(close_contacts): " + str(len(close_contacts))
      print "contacts"
      for key, val in contacts.iteritems():
        if val != 0:
          print (key, val)
      print "len(contacts): " + str(len(contacts))
      # The keys of these dicts are pairs of atomtypes, but the keys are
      # sorted so that ("C", "O") is always written as "C_O". Thus, for N
      # atom types, there are N*(N+1)/2 unique pairs.
      num_atoms = len(Binana.atom_types)
      assert len(close_contacts) == num_atoms*(num_atoms+1)/2
      assert len(contacts) == num_atoms*(num_atoms+1)/2

  def test_compute_pi_pi_stacking(self):
    """
    TestBinana: Compute Pi-Pi Stacking.
    """
    # 1zea is the only example that has any pi-stacking.
    pi_stacking_dict = {}
    for name, protein, ligand in self.test_cases:
      pi_stacking_dict[name] = compute_pi_pi_stacking(
          ligand, protein)
    for name, pi_stacking in pi_stacking_dict.iteritems():
      print "Processing pi-stacking for %s" % name
      assert len(pi_stacking) == 3
      print pi_stacking
      assert "STACKING_ALPHA" in pi_stacking
      assert "STACKING_BETA" in pi_stacking
      assert "STACKING_OTHER" in pi_stacking


  def test_compute_pi_t(self):
    """
    TestBinana: Compute Pi-T Interactions.

    TODO(rbharath): I believe that the imatininb-cabl complex has a pi-T
    interaction. This code has a bug since it reports that no such
    interaction is found.
    """
    pi_t_dict = {}
    for name, protein, ligand in self.test_cases:
      pi_t_dict[name] = compute_pi_t(
          ligand, protein)
    for name, pi_t in pi_t_dict.iteritems():
      print "Processing pi-T for %s" % name
      assert len(pi_t) == 3
      assert "T-SHAPED_ALPHA" in pi_t
      assert "T-SHAPED_BETA" in pi_t
      assert "T-SHAPED_OTHER" in pi_t

  def test_compute_pi_cation(self):
    """
    TestBinana: Compute Pi-Cation Interactions.
    """
    pi_cation_dict = {}
    for name, protein, ligand in self.test_cases:
      pi_cation_dict[name] = compute_pi_cation(
          ligand, protein)
    for name, pi_cation in pi_cation_dict.iteritems():
      print "Processing pi-cation for %s" % name
      assert len(pi_cation) == 6
      assert 'PI-CATION_LIGAND-CHARGED_ALPHA' in pi_cation
      assert 'PI-CATION_LIGAND-CHARGED_BETA' in pi_cation
      assert 'PI-CATION_LIGAND-CHARGED_OTHER' in pi_cation
      assert 'PI-CATION_RECEPTOR-CHARGED_ALPHA' in pi_cation
      assert 'PI-CATION_RECEPTOR-CHARGED_BETA' in pi_cation
      assert 'PI-CATION_RECEPTOR-CHARGED_OTHER' in pi_cation

  def test_compute_salt_bridges(self):
    """
    TestBinana: Compute Salt Bridges.

    TODO(bramsundar): None of the examples contain salt-bridge interactions. Find a
    complex with an actual salt-bridge interaction.
    """
    salt_bridges_dict = {}
    for name, protein, ligand in self.test_cases:
      salt_bridges_dict[name] = compute_salt_bridges(
          ligand, protein)
    for name, salt_bridges in salt_bridges_dict.iteritems():
      print "Processing salt-bridges for %s" % name
      assert len(salt_bridges) == 3
      print salt_bridges
      assert 'SALT-BRIDGE_ALPHA' in salt_bridges
      assert 'SALT-BRIDGE_BETA' in salt_bridges
      assert 'SALT-BRIDGE_OTHER' in salt_bridges

  def test_compute_input_vector(self):
    """
    TestBinana: Compute Input Vector.
    """
    features_dict = {}
    for name, protein, ligand in self.test_cases:
      features_dict[name] = self.binana.compute_input_vector(
          ligand, protein)
    num_atoms = len(Binana.atom_types)
    # Lengths:
    # ligand_receptor_close_contacts: N*(N+1)/2
    # ligand_receptor_contacts: N*(N+1)/2
    # ligand_receptor_electrostatics: N*(N+1)/2
    # ligand_atom_counts: N
    # hbonds: 12
    # hydrophobics: 6
    # stacking: 3
    # pi_cation: 6
    # t_shaped: 3
    # active_site_flexibility: 6
    # salt_bridges: 3
    # rotatable_boonds_count: 1
    total_len = (3*num_atoms*(num_atoms+1)/2 + num_atoms + 12 + 6 + 3 + 6 +
        3 + 6 + 3 + 1)
    for name, input_vector in features_dict.iteritems():
      print "Processing input-vector for %s" % name
      assert len(input_vector) == total_len
