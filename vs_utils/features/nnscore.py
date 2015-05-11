"""
The following code implements a featurizer based on NNScore 2.0.1

## The following notice is copied from the original NNScore file.
# NNScore 2.01 is released under the GNU General Public License (see
# http://www.gnu.org/licenses/gpl.html).
# If you have any questions, comments, or suggestions, please don't
# hesitate to contact me, Jacob Durrant, at jdurrant [at] ucsd [dot]
# edu. If you use NNScore 2.01 in your work, please cite [REFERENCE
# HERE].
"""

__author__ = "Bharath Ramsundar"
__license__ = "GNU General Public License"

import textwrap
import math
import os
import sys
import textwrap
import glob
import cPickle
import numpy as np
from vs_utils.features import Featurizer
from vs_utils.features.nnscore_helper import MathFunctions
from vs_utils.features.nnscore_helper import PDB
from vs_utils.features.nnscore_helper import Point

# TODO(bramsundar): Many places in this file use long switch
# statements. Could there be a cleaner way to structure this file?

# TODO(bramsundar): Add some tests for the classes and functions in
# this file. Will help verify that cleanup didn't bork anything.

# TODO(bramsundar): How does vs_utils handle scripts? Should I factor
# the script part of this file elsewhere?

# TODO(bramsundar): Many arbitrary choices in cutoff-values. I think
# that moving to a 3D-convnet might actually be necessary to remove
# the high-degree of hand-tuning here.

ELECTROSTATIC_JOULE_PER_MOL = 138.94238460104697e4 # units?
# O-H distance is 0.96 A, N-H is 1.01 A. See
# http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html
H_BOND_DIST = 1.3 # angstroms
H_BOND_ANGLE = 40 # degrees
# If atoms are < 2.5 A apart, we count it as a close contact
CLOSE_CONTACT_CUTOFF = 2.5
# If receptor and ligand atoms are > 4 A apart, we consider them
# unable to interact with simple electrostatics.
CONTACT_CUTOFF = 4 # angstroms
# "PI-Stacking Interactions ALIVE AND WELL IN PROTEINS" says
# distance of 7.5 A is good cutoff. This seems really big to me,
# except that pi-pi interactions (parallel) are actually usually
# off centered. Interesting paper.  Note that adenine and
# tryptophan count as two aromatic rings. So, for example, an
# interaction between these two, if positioned correctly, could
# count for 4 pi-pi interactions.
PI_PI_CUTOFF = 7.5
# Cation-pi interaction cutoff based on
# "Cation-pi interactions in structural biology."
CATION_PI_CUTOFF = 6.0
# 4  is good cutoff for salt bridges according to
# "Close-Range Electrostatic Interactions in Proteins",
# but looking at complexes, I decided to go with 5.5 A
SALT_BRIDGE_CUTOFF = 5.5
# This is perhaps controversial. I noticed that often a pi-cation
# interaction or other pi interaction was only slightly off, but
# looking at the structure, it was clearly supposed to be a pi-cation
# interaction. I've decided then to artificially expand the radius of
# each pi ring. Think of this as adding in a VDW radius, or
# accounting for poor crystal-structure resolution, or whatever you
# want to justify it.
PI_PADDING = 0.75


def hashtable_entry_add_one(hashtable, key, toadd = 1):
  # note that dictionaries (hashtables) are passed by reference in python
  if hashtable.has_key(key):
    hashtable[key] = hashtable[key] + toadd
  else:
    hashtable[key] = toadd


def extend_list_by_dictionaries(lst, dictionaries):
  newlist = lst[:]
  for dictionary in dictionaries:
    # first, sort the dictionary by the key
    keys = dictionary.keys()
    keys.sort()

    # now make a list of the values
    vals = []
    for key in keys:
      vals.append(dictionary[key])

    # now append vals to the list
    newlist.extend(vals)

  # return the extended list
  return newlist

def center(string, length):
  while len(string) < length:
    string = " " + string
    if len(string) < length:
      string = string + " "
  return string


class binana:
  """
  binana extracts a fingerprint from a provided binding pose.

  The binana feature vector transforms a ligand-receptor binding pose
  into a feature vector. The feature vector has the following
  components:

    -vina_output: Components of vina's score function.
    -ligand_receptor_contacts: List of contacts between ligand and
       receptor atoms (< 4 A)
    -ligand_receptor_electrostatics: Coulomb energy between contacting
       ligand and receptor atoms.
    -ligand_atom_types: The atom types in the ligand.
    -ligand_receptor_close_contacts: List of close contacts between
       ligand and receptor (< 2.5 A)
    -hbonds: List of hydrogen bonds.
    -hydrophobic: List of hydrophobic contacts.
    -stacking: List of pi-pi stacking.
    -pi_cation: List of pi-cation interactions.
    -t_shaped: List of T-shaped interactions.
      TODO(bramsundar): What are T-shaped interactions?
    -active_site_flexibility: Considers whether the receptor atoms are
       backbone or sidechain and whether they are part of
       alpha-helices or beta-sheets.
    -salt_bridges: List of salt-bridges between ligand and receptor.
    -rotatable_bonds_count: Count of (ligand(?), receptor(?))
       rotatable bonds.
  """

  # supporting functions
  def get_vina_output(self):
    """Invoke vina on ligand and receptor."""
    # Now save the files
    preface ="REMARK "

    # Now get vina
    vina_output = getCommandOutput2(parameters.params['vina_executable'] +
        ' --score_only --receptor ' + receptor_pdbqt_filename + ' --ligand '
        + ligand_pdbqt_filename)

    #print vina_output
    vina_output = vina_output.split("\n")
    vina_affinity = 0.0
    vina_gauss_1 = 0.0
    vina_gauss_2 = 0.0
    vina_repulsion = 0.0
    vina_hydrophobic = 0.0
    vina_hydrogen = 0.0
    for item in vina_output:
      item = item.strip()
      if "Affinity" in item:
        vina_affinity = float(item.replace("Affinity: ","").replace(" (kcal/mol)",""))
      if "gauss 1" in item:
        vina_gauss_1 = float(item.replace("gauss 1     : ",""))
      if "gauss 2" in item:
        vina_gauss_2 = float(item.replace("gauss 2     : ",""))
      if "repulsion" in item:
        vina_repulsion = float(item.replace("repulsion   : ",""))
      if "hydrophobic" in item:
        vina_hydrophobic = float(item.replace("hydrophobic : ",""))
      if "Hydrogen" in item:
        vina_hydrogen = float(item.replace("Hydrogen    : ",""))

    vina_output = [vina_affinity, vina_gauss_1, vina_gauss_2,
        vina_repulsion, vina_hydrophobic, vina_hydrogen]

  def compute_hydrophobic_contacts(self, ligand, receptor):
    """
    Compute possible hydrophobic contacts between ligand and atom.

    Returns a dictionary whose keys are atompairs of type
    "${RESIDUETYPE}_${RECEPTOR_ATOM}" where RESIDUETYPE is either "SIDECHAIN" or
    "BACKBONE" and RECEPTOR_ATOM is "O" or "C" or etc. The
    values count the number of hydrophobic contacts.

    Parameters
    ----------
    ligand: PDB
      A PDB Object describing the ligand molecule.
    receptor: PDB
      A PDB object describing the receptor protein.

    """
    # Now see if there's hydrophobic contacts (C-C contacts)
    hydrophobics = {
      'BACKBONE_ALPHA': 0, 'BACKBONE_BETA': 0, 'BACKBONE_OTHER': 0,
      'SIDECHAIN_ALPHA': 0, 'SIDECHAIN_BETA': 0, 'SIDECHAIN_OTHER': 0
      }
    for ligand_atom_index in ligand.all_atoms:
      ligand_atom = ligand.all_atoms[ligand_atom_index]
      for receptor_atom_index in receptor.all_atoms:
        receptor_atom = receptor.all_atoms[receptor_atom_index]
        if ligand_atom.element == "C" and receptor_atom.element == "C":
          hydrophobic_key = (receptor_atom.side_chain_or_backbone() +
              "_" + receptor_atom.structure)

          hashtable_entry_add_one(hydrophobics, hydrophobic_key)
    return hydrophobics

  def compute_electrostatic_energy(self, ligand, receptor):
    """
    Compute electrostatic energy between ligand and atom.

    Returns a dictionary whose keys are atompairs of type
    "${RECEPTOR_ATOM}_${LIGAND_ATOM}". The ATOM terms can equal "C", "O",
    etc. The values are the associated coulomb energies for this pair.

    Parameters
    ----------
    ligand: PDB
      A PDB Object describing the ligand molecule.
    receptor: PDB
      A PDB object describing the receptor protein.

    """
    ligand_receptor_electrostatics = {
        "A_MG": 0, "A_MN": 0, "BR_SA": 0, "CL_FE": 0, "CL_MG": 0,
        "CL_MN": 0, "CL_NA": 0, "CL_P": 0, "CL_S": 0, "CL_ZN": 0,
        "CU_HD": 0, "CU_N": 0, "FE_NA": 0, "FE_SA": 0, "MG_N": 0,
        "MG_S": 0, "MG_SA": 0, "MN_NA": 0, "MN_S": 0, "MN_SA": 0,
        "NA_P": 0, "P_S": 0, "P_SA": 0, "S_SA": 0, "A_A": 0.0, "A_BR":
        0.0, "A_C": 0.0, "A_CL": 0.0, "A_F": 0.0, "A_FE": 0.0, "A_HD":
        0.0, "A_I": 0.0, "A_N": 0.0, "A_NA": 0.0, "A_OA": 0.0, "A_P":
        0.0, "A_S": 0.0, "A_SA": 0.0, "A_ZN": 0.0, "BR_C": 0.0,
        "BR_HD": 0.0, "BR_N": 0.0, "BR_OA": 0.0, "C_C": 0.0, "C_CL":
        0.0, "C_F": 0.0, "C_FE": 0.0, "C_HD": 0.0, "C_I": 0.0,
        "CL_HD": 0.0, "CL_N": 0.0, "CL_OA": 0.0, "CL_SA": 0.0, "C_MG":
        0.0, "C_MN": 0.0, "C_N": 0.0, "C_NA": 0.0, "C_OA": 0.0, "C_P":
        0.0, "C_S": 0.0, "C_SA": 0.0, "C_ZN": 0.0, "FE_HD": 0.0,
        "FE_N": 0.0, "FE_OA": 0.0, "F_HD": 0.0, "F_N": 0.0, "F_OA":
        0.0, "F_SA": 0.0, "HD_HD": 0.0, "HD_I": 0.0, "HD_MG": 0.0,
        "HD_MN": 0.0, "HD_N": 0.0, "HD_NA": 0.0, "HD_OA": 0.0, "HD_P":
        0.0, "HD_S": 0.0, "HD_SA": 0.0, "HD_ZN": 0.0, "I_N": 0.0,
        "I_OA": 0.0, "MG_NA": 0.0, "MG_OA": 0.0, "MG_P": 0.0, "MN_N":
        0.0, "MN_OA": 0.0, "MN_P": 0.0, "NA_OA": 0.0, "NA_S": 0.0,
        "NA_SA": 0.0, "NA_ZN": 0.0, "N_N": 0.0, "N_NA": 0.0, "N_OA":
        0.0, "N_P": 0.0, "N_S": 0.0, "N_SA": 0.0, "N_ZN": 0.0,
        "OA_OA": 0.0, "OA_P": 0.0, "OA_S": 0.0, "OA_SA": 0.0, "OA_ZN":
        0.0, "P_ZN": 0.0, "SA_SA": 0.0, "SA_ZN": 0.0, "S_ZN": 0,
        "F_ZN": 0}
    for ligand_atom_index in ligand.all_atoms:
      ligand_atom = ligand.all_atoms[ligand_atom_index]
      for receptor_atom_index in receptor.all_atoms:
        receptor_atom = receptor.all_atoms[receptor_atom_index]
        atomtypes = [ligand_atom.atomtype, receptor_atom.atomtype]
        atomstr = "_".join(sorted(atomtypes))
        dist = ligand_atom.coordinates.dist_to(receptor_atom.coordinates)
        # calculate electrostatic energies for all less than 4 A
        ligand_charge = ligand_atom.charge
        receptor_charge = receptor_atom.charge
        # to convert into J/mol; might be nice to double check this
        # TODO(bramsundar): What are units of
        # ligand_charge/receptor_charge?
        coulomb_energy = ((ligand_charge * receptor_charge / dist)
            * ELECTROSTATIC_JOULE_PER_MOL)
        hashtable_entry_add_one(ligand_receptor_electrostatics,
            atomstr, coulomb_energy)
    return ligand_receptor_electrostatics

  def compute_active_site_flexibility(self, ligand, receptor):
    """
    Compute statistics to judge active-site flexibility

    Returns a dictionary whose keys are atompairs of type
    "${RESIDUETYPE}_${STRUCTURE}" where RESIDUETYPE is either "SIDECHAIN" or
    "BACKBONE" and STRUCTURE is either ALPHA, BETA, or OTHER and
    corresponds to the protein secondary structure of the current residue.

    Parameters
    ----------
    ligand: PDB
      A PDB Object describing the ligand molecule.
    receptor: PDB
      A PDB object describing the receptor protein.

    """
    active_site_flexibility = {
      'BACKBONE_ALPHA': 0, 'BACKBONE_BETA': 0, 'BACKBONE_OTHER': 0,
      'SIDECHAIN_ALPHA': 0, 'SIDECHAIN_BETA': 0, 'SIDECHAIN_OTHER': 0
      }
    for receptor_atom_index in receptor.all_atoms:
      receptor_atom = receptor.all_atoms[receptor_atom_index]

      # TODO(rbharath): This feels like a silent failure mode... Introduce
      # an upstream fix.
      if receptor_atom.structure == "":
        structure = "OTHER"
      else:
        structure = receptor_atom.structure
      flexibility_key = (receptor_atom.side_chain_or_backbone() + "_"
          + structure)
      hashtable_entry_add_one(active_site_flexibility, flexibility_key)
    return active_site_flexibility


  def compute_hydrogen_bonds(self, ligand, receptor):
    """
    Computes hydrogen bonds between ligand and receptor.

    Returns a dictionary whose keys are of form
    HDONOR-${COMMENT}_${RESIDUETYPE}_${STRUCTURE} where COMMENT is either
    "RECEPTOR" or "LIGAND", RESIDUETYPE is "BACKBONE" or "SIDECHAIN" and
    where STRUCTURE is "ALPHA" or "BETA" or "OTHER". The values are counts
    of the numbers of hydrogen bonds associated with the given keys.

    Parameters
    ----------
    ligand: PDB
      A PDB Object describing the ligand molecule.
    receptor: PDB
      A PDB object describing the receptor protein.
    """
    hbonds = {
      'HDONOR-LIGAND_BACKBONE_ALPHA': 0,
      'HDONOR-LIGAND_BACKBONE_BETA': 0,
      'HDONOR-LIGAND_BACKBONE_OTHER': 0,
      'HDONOR-LIGAND_SIDECHAIN_ALPHA': 0,
      'HDONOR-LIGAND_SIDECHAIN_BETA': 0,
      'HDONOR-LIGAND_SIDECHAIN_OTHER': 0,
      'HDONOR-RECEPTOR_BACKBONE_ALPHA': 0,
      'HDONOR-RECEPTOR_BACKBONE_BETA': 0,
      'HDONOR-RECEPTOR_BACKBONE_OTHER': 0,
      'HDONOR-RECEPTOR_SIDECHAIN_ALPHA': 0,
      'HDONOR-RECEPTOR_SIDECHAIN_BETA': 0,
      'HDONOR-RECEPTOR_SIDECHAIN_OTHER': 0}
    for ligand_atom_index in ligand.all_atoms:
      ligand_atom = ligand.all_atoms[ligand_atom_index]
      for receptor_atom_index in receptor.all_atoms:
        receptor_atom = receptor.all_atoms[receptor_atom_index]
        # Now see if there's some sort of hydrogen bond between
        # these two atoms. distance cutoff = H_BOND_DIST, angle cutoff =
        # H_BOND_ANGLE.
        # Note that this is liberal.
        if ((ligand_atom.element == "O" or ligand_atom.element == "N")
          and (receptor_atom.element == "O"
              or receptor_atom.element == "N")):

          # now build a list of all the hydrogens close to these
          # atoms
          hydrogens = []

          for atm_index in ligand.all_atoms:
            if ligand.all_atoms[atm_index].element == "H":
              # so it's a hydrogen
              if (ligand.all_atoms[atm_index].coordinates.dist_to(
                  ligand_atom.coordinates) < H_BOND_DIST):
                ligand.all_atoms[atm_index].comment = "LIGAND"
                hydrogens.append(ligand.all_atoms[atm_index])

          for atm_index in receptor.all_atoms:
            if receptor.all_atoms[atm_index].element == "H": # so it's a hydrogen
              if (receptor.all_atoms[atm_index].coordinates.dist_to(
                  receptor_atom.coordinates) < H_BOND_DIST):
                receptor.all_atoms[atm_index].comment = "RECEPTOR"
                hydrogens.append(receptor.all_atoms[atm_index])

          # now we need to check the angles
          # TODO(bramsundar): Verify this heuristic and add a source
          # citation.
          for hydrogen in hydrogens:
            if (math.fabs(180 - self.functions.angle_between_three_points(
                  ligand_atom.coordinates,
                  hydrogen.coordinates,
                  receptor_atom.coordinates) * 180.0 / math.pi) <=
                  H_BOND_ANGLE):
              # TODO(rbharath): This feels like a silent failure mode... Introduce
              # an upstream fix.
              if receptor_atom.structure == "":
                structure = "OTHER"
              else:
                structure = receptor_atom.structure
              hbonds_key = ("HDONOR-" + hydrogen.comment + "_" +
                  receptor_atom.side_chain_or_backbone() + "_" +
                  structure)
              hashtable_entry_add_one(hbonds, hbonds_key)
    return hbonds

  def compute_ligand_atom_counts(self, ligand):
    """Counts atoms of each type in given ligand.

    Returns a dictionary that maps atom types ("C", "O", etc.) to
    counts.

    Parameters
    ----------
    ligand: PDB Object
      Stores ligand information.

    Returns
    -------
    ligand_atom_types: dictionary
      Keys are atom types; values are integer counts.
    """
    # TODO(rbharath): What is atom type A here?
    ligand_atom_types = {
        "A": 0, "BR": 0, "C": 0, "CL": 0, "F": 0, "H": 0,  "HD": 0, "I": 0,
        "N": 0, "NA": 0, "O": 0, "OA": 0, "P": 0, "S": 0, "SA": 0}
    for ligand_atom_index in ligand.all_atoms:
      ligand_atom = ligand.all_atoms[ligand_atom_index]
      hashtable_entry_add_one(ligand_atom_types, ligand_atom.atomtype)
    return ligand_atom_types

  def compute_ligand_receptor_contacts(self, ligand, receptor):
    """Compute distance measurements for ligand-receptor atom pairs.

    Returns two dictionaries, each of whose keys are of form
    ATOMTYPE_ATOMTYPE.

    Parameters
    ----------
    ligand: TODO(bramsundar)
    receptor: TODO(bramsundar)
    """
    ligand_receptor_close_contacts = {
        "A_A": 0, "A_C": 0, "A_CL": 0, "A_F": 0, "A_FE": 0, "A_MG": 0,
        "A_MN": 0, "A_NA": 0, "A_SA": 0, "BR_C": 0, "BR_OA": 0, "C_CL":
        0, "CD_OA": 0, "CL_FE": 0, "CL_MG": 0, "CL_N": 0, "CL_OA": 0,
        "CL_ZN": 0, "C_MN": 0, "C_NA": 0, "F_N": 0, "F_SA": 0, "F_ZN":
        0, "HD_MN": 0, "MN_N": 0, "NA_SA": 0, "N_SA": 0, "A_HD": 0,
        "A_N": 0, "A_OA": 0, "A_ZN": 0, "BR_HD": 0, "C_C": 0, "C_F": 0,
        "C_HD": 0, "CL_HD": 0, "C_MG": 0, "C_N": 0, "C_OA": 0, "C_SA":
        0, "C_ZN": 0, "FE_HD": 0, "FE_N": 0, "FE_OA": 0, "F_HD": 0,
        "F_OA": 0, "HD_HD": 0, "HD_I": 0, "HD_MG": 0, "HD_N": 0,
        "HD_NA": 0, "HD_OA": 0, "HD_P": 0, "HD_S": 0, "HD_SA": 0,
        "HD_ZN": 0, "MG_NA": 0, "MG_OA": 0, "MN_OA": 0, "NA_OA": 0,
        "NA_ZN": 0, "N_N": 0, "N_NA": 0, "N_OA": 0, "N_ZN": 0, "OA_OA":
        0, "OA_SA": 0, "OA_ZN": 0, "SA_ZN": 0, "S_ZN": 0}
    ligand_receptor_contacts = {
        "A_CU": 0, "A_MG": 0, "A_MN": 0, "BR_SA": 0, "C_CD": 0,
        "CL_FE": 0, "CL_MG": 0, "CL_MN": 0, "CL_NA": 0, "CL_P": 0,
        "CL_S": 0, "CL_ZN": 0, "CU_HD": 0, "CU_N": 0, "FE_NA": 0,
        "FE_SA": 0, "MG_N": 0, "MG_S": 0, "MG_SA": 0, "MN_NA": 0,
        "MN_S": 0, "MN_SA": 0, "NA_P": 0, "P_S": 0, "P_SA": 0, "S_SA":
        0, "A_A": 0, "A_BR": 0, "A_C": 0, "A_CL": 0, "A_F": 0, "A_FE":
        0, "A_HD": 0, "A_I": 0, "A_N": 0, "A_NA": 0, "A_OA": 0, "A_P":
        0, "A_S": 0, "A_SA": 0, "A_ZN": 0, "BR_C": 0, "BR_HD": 0,
        "BR_N": 0, "BR_OA": 0, "C_C": 0, "C_CL": 0, "C_F": 0, "C_FE":
        0, "C_HD": 0, "C_I": 0, "CL_HD": 0, "CL_N": 0, "CL_OA": 0,
        "CL_SA": 0, "C_MG": 0, "C_MN": 0, "C_N": 0, "C_NA": 0, "C_OA":
        0, "C_P": 0, "C_S": 0, "C_SA": 0, "C_ZN": 0, "FE_HD": 0,
        "FE_N": 0, "FE_OA": 0, "F_HD": 0, "F_N": 0, "F_OA": 0, "F_SA":
        0, "HD_HD": 0, "HD_I": 0, "HD_MG": 0, "HD_MN": 0, "HD_N": 0,
        "HD_NA": 0, "HD_OA": 0, "HD_P": 0, "HD_S": 0, "HD_SA": 0,
        "HD_ZN": 0, "I_N": 0, "I_OA": 0, "MG_NA": 0, "MG_OA": 0,
        "MG_P": 0, "MN_N": 0, "MN_OA": 0, "MN_P": 0, "NA_OA": 0,
        "NA_S": 0, "NA_SA": 0, "NA_ZN": 0, "N_N": 0, "N_NA": 0,
        "N_OA": 0, "N_P": 0, "N_S": 0, "N_SA": 0, "N_ZN": 0, "OA_OA":
        0, "OA_P": 0, "OA_S": 0, "OA_SA": 0, "OA_ZN": 0, "P_ZN": 0,
        "SA_SA": 0, "SA_ZN": 0, "S_ZN": 0}
    for ligand_atom_index in ligand.all_atoms:
      for receptor_atom_index in receptor.all_atoms:
        ligand_atom = ligand.all_atoms[ligand_atom_index]
        receptor_atom = receptor.all_atoms[receptor_atom_index]

        dist = ligand_atom.coordinates.dist_to(receptor_atom.coordinates)
        atomtypes = [ligand_atom.atomtype, receptor_atom.atomtype]
        atomstr = "_".join(sorted(atomtypes))
        if dist < CONTACT_CUTOFF:
          hashtable_entry_add_one(ligand_receptor_contacts,
              atomstr)
        if dist < CLOSE_CONTACT_CUTOFF:
          hashtable_entry_add_one(ligand_receptor_close_contacts,
              atomstr)
    return ligand_receptor_close_contacts, ligand_receptor_contacts

  def compute_pi_pi_stacking(self, ligand, receptor):
    """
    Computes pi-pi interactions.

    Parameters
    ----------
    ligand: PDB Object.
      small molecule to dock.
    receptor: PDB Object
      protein to dock agains.
    """
    pi_stacking = {'ALPHA': 0, 'BETA': 0, 'OTHER': 0}
    for lig_aromatic in ligand.aromatic_rings:
      for rec_aromatic in receptor.aromatic_rings:
        dist = lig_aromatic.center.dist_to(rec_aromatic.center)
        if dist < PI_PI_CUTOFF:
          # so there could be some pi-pi interactions.  Now, let's
          # check for stacking interactions. Are the two pi's roughly
          # parallel?
          lig_aromatic_norm_vector = Point(coords=np.array([lig_aromatic.plane_coeff[0],
              lig_aromatic.plane_coeff[1], lig_aromatic.plane_coeff[2]]))
          rec_aromatic_norm_vector = Point(coords=np.array([rec_aromatic.plane_coeff[0],
              rec_aromatic.plane_coeff[1], rec_aromatic.plane_coeff[2]]))
          angle_between_planes = self.functions.angle_between_points(
              lig_aromatic_norm_vector, rec_aromatic_norm_vector) * 180.0/math.pi

          if math.fabs(angle_between_planes-0) < 30.0 or math.fabs(angle_between_planes-180) < 30.0:
            # so they're more or less parallel, it's probably pi-pi
            # stackingoutput_dir now, pi-pi are not usually right on
            # top of each other. They're often staggered. So I don't
            # want to just look at the centers of the rings and
            # compare. Let's look at each of the atoms.  do atom of
            # the atoms of one ring, when projected onto the plane of
            # the other, fall within that other ring?

            # start by assuming it's not a pi-pi stacking interaction
            pi_pi = False
            for ligand_ring_index in lig_aromatic.indices:
              # project the ligand atom onto the plane of the receptor ring
              pt_on_receptor_plane = self.functions.project_point_onto_plane(
                  ligand.all_atoms[ligand_ring_index].coordinates,
                  rec_aromatic.plane_coeff)
              if (pt_on_receptor_plane.dist_to(rec_aromatic.center)
                 <= rec_aromatic.radius + PI_PADDING):
                pi_pi = True
                break

            if pi_pi == False:
              # if you've already determined it's a pi-pi stacking interaction, no need to keep trying
              for receptor_ring_index in rec_aromatic.indices:
                # project the ligand atom onto the plane of the receptor ring
                pt_on_ligand_plane = self.functions.project_point_onto_plane(
                    receptor.all_atoms[receptor_ring_index].coordinates,
                    lig_aromatic.plane_coeff)
                if (pt_on_ligand_plane.dist_to(lig_aromatic.center)
                    <= lig_aromatic.radius + PI_PADDING):
                  pi_pi = True
                  break

            if pi_pi == True:
                structure = receptor.all_atoms[rec_aromatic.indices[0]].structure
                if structure == "":
                  # since it could be interacting with a cofactor or something
                  structure = "OTHER"
                key = "STACKING_" + structure


                hashtable_entry_add_one(pi_stacking, key)
    return pi_stacking

  def compute_pi_T(self, ligand, receptor):
    """
    Computes T-shaped pi-pi interactions.

    Parameters
    ----------
    ligand: PDB Object.
      small molecule to dock.
    receptor: PDB Object
      protein to dock agains.
    """
    pi_T = {'ALPHA': 0, 'BETA': 0, 'OTHER': 0}
    for lig_aromatic in ligand.aromatic_rings:
      for rec_aromatic in receptor.aromatic_rings:
        lig_aromatic_norm_vector = Point(coords=np.array([lig_aromatic.plane_coeff[0],
            lig_aromatic.plane_coeff[1], lig_aromatic.plane_coeff[2]]))
        rec_aromatic_norm_vector = Point(coords=np.array([rec_aromatic.plane_coeff[0],
            rec_aromatic.plane_coeff[1], rec_aromatic.plane_coeff[2]]))
        angle_between_planes = self.functions.angle_between_points(
            lig_aromatic_norm_vector, rec_aromatic_norm_vector) * 180.0/math.pi
        if math.fabs(angle_between_planes-90) < 30.0 or math.fabs(angle_between_planes-270) < 30.0:
          # so they're more or less perpendicular, it's probably a
          # pi-edge interaction having looked at many structures, I
          # noticed the algorithm was identifying T-pi reactions
          # when the two rings were in fact quite distant, often
          # with other atoms in between. Eye-balling it, requiring
          # that at their closest they be at least 5 A apart seems
          # to separate the good T's from the bad
          min_dist = 100.0
          for ligand_ind in lig_aromatic.indices:
            ligand_at = ligand.all_atoms[ligand_ind]
            for receptor_ind in rec_aromatic.indices:
              receptor_at = receptor.all_atoms[receptor_ind]
              dist = ligand_at.coordinates.dist_to(receptor_at.coordinates)
              if dist < min_dist: min_dist = dist

          if min_dist <= 5.0:
            # so at their closest points, the two rings come within
            # 5 A of each other.

            # okay, is the ligand pi pointing into the receptor
            # pi, or the other way around?  first, project the
            # center of the ligand pi onto the plane of the
            # receptor pi, and vs. versa

            # This could be directional somehow, like a hydrogen
            # bond.

            pt_on_receptor_plane = self.functions.project_point_onto_plane(
                lig_aromatic.center, rec_aromatic.plane_coeff)
            pt_on_ligand_plane = self.functions.project_point_onto_plane(
                rec_aromatic.center, lig_aromatic.plane_coeff)

            # now, if it's a true pi-T interaction, this projected
            # point should fall within the ring whose plane it's
            # been projected into.
            if ((pt_on_receptor_plane.dist_to(rec_aromatic.center)
              <= rec_aromatic.radius + PI_PADDING)
              or (pt_on_ligand_plane.dist_to(lig_aromatic.center)
              <= lig_aromatic.radius + PI_PADDING)):

              # so it is in the ring on the projected plane.
              structure = receptor.all_atoms[rec_aromatic.indices[0]].structure
              if structure == "":
                # since it could be interacting with a cofactor or something
                structure = "OTHER"
              key = "T-SHAPED_" + structure

              hashtable_entry_add_one(pi_T, key)
    return pi_T

  def compute_pi_cation(self, ligand, receptor):
    """
    Computes number of pi-cation interactions.

    Returns a dictionary whose keys are of form
    ${MOLTYPE}-CHARGED_${STRUCTURE} where MOLTYPE is either "LIGAND" or
    "RECEPTOR" and STRUCTURE is "ALPHA" or "BETA" or "OTHER".

    Parameters
    ----------
    ligand: PDB Object
      small molecule to dock.
    receptor: PDB Object
      protein to dock agains.
    """
    pi_cation = {
      'LIGAND-CHARGED_ALPHA': 0, 'LIGAND-CHARGED_BETA': 0,
      'LIGAND-CHARGED_OTHER': 0, 'RECEPTOR-CHARGED_ALPHA': 0,
      'RECEPTOR-CHARGED_BETA': 0, 'RECEPTOR-CHARGED_OTHER': 0}
    for aromatic in receptor.aromatic_rings:
      for charged in ligand.charges:
        if charged.positive == True: # so only consider positive charges
          if (charged.coordinates.dist_to(aromatic.center) < CATION_PI_CUTOFF):

            # project the charged onto the plane of the aromatic
            charge_projected = self.functions.project_point_onto_plane(
                charged.coordinates,aromatic.plane_coeff)
            if (charge_projected.dist_to(aromatic.center)
              < aromatic.radius + PI_PADDING):
              structure = receptor.all_atoms[aromatic.indices[0]].structure
              if structure == "":
                # since it could be interacting with a cofactor or something
                structure = "OTHER"
              key = "PI-CATION_LIGAND-CHARGED_" + structure

              hashtable_entry_add_one(pi_cation, key)

    for aromatic in ligand.aromatic_rings:
      # now it's the ligand that has the aromatic group
      for charged in receptor.charges:
        if charged.positive == True: # so only consider positive charges
          if (charged.coordinates.dist_to(aromatic.center) < CATION_PI_CUTOFF):
            charge_projected = self.functions.project_point_onto_plane(
                charged.coordinates,aromatic.plane_coeff)
            if charge_projected.dist_to(aromatic.center) < aromatic.radius + PI_PADDING:
              structure = receptor.all_atoms[charged.indices[0]].structure
              if structure == "":
                # since it could be interacting with a cofactor or something
                structure = "OTHER"
              key = "PI-CATION_RECEPTOR-CHARGED_" + structure

              hashtable_entry_add_one(pi_cation, key)
    return pi_cation

  def compute_salt_bridges(self, ligand, receptor):
    """
    Computes number of ligand-receptor salt bridges.

    Returns a dictionary with keys of form SALT-BRIDGE_${STRUCTURE} where
    STRUCTURE is "ALPHA" or "BETA" or "OTHER."

    Parameters
    ----------
    ligand: PDB Object
      small molecule to dock.
    receptor: PDB Object
      protein to dock agains.
    """
    salt_bridges = {'ALPHA': 0, 'BETA': 0, 'OTHER': 0}
    for receptor_charge in receptor.charges:
      for ligand_charge in ligand.charges:
        if ligand_charge.positive != receptor_charge.positive:
          # so they have oppositve charges
          if (ligand_charge.coordinates.dist_to(
              receptor_charge.coordinates) < SALT_BRIDGE_CUTOFF):
            structure = receptor.all_atoms[receptor_charge.indices[0]].structure
            if structure == "":
              # since it could be interacting with a cofactor or something
              structure = "OTHER"
            key = "SALT-BRIDGE_" + structure

            hashtable_entry_add_one(salt_bridges, key)
    return salt_bridges

  def compute_input_vector_from_files(self, ligand_pdbqt_filename,
      receptor_pdbqt_filename, line_header):
    """Computes feature vector for ligand-receptor pair.

    Parameters
    ----------
    ligand_pdbqt_filename: string
      path to ligand's pdbqt file.
    receptor_pdbqt_filename: string
      path to receptor pdbqt file.
    line_header: string
      TODO(bramsundar): line separator in PDB files (?)
    """
    # Load receptor and ligand from file.
    receptor = PDB()
    receptor.load_PDB_from_file(receptor_pdbqt_filename, line_header)
    receptor.assign_secondary_structure()
    ligand = PDB()
    ligand.load_PDB_from_file(ligand_pdbqt_filename, line_header)
    self.compute_input_vector(ligand, pdb)

  def compute_input_vector(self, ligand, receptor):
    """Computes feature vector for ligand-receptor pair.

    Parameters
    ----------
    ligand: PDB object
      Contains loaded ligand.
    receptor: PDB object
      Contains loaded receptor.
    """

    rotatable_bonds_count = {'rot_bonds': ligand.rotatable_bonds_count}
    ligand_receptor_close_contacts, ligand_receptor_contacts = (
      self.compute_ligand_receptor_contacts(ligand, receptor))
    ligand_receptor_electrostatics = (
      self.compute_electrostatic_energy(ligand, receptor))
    ligand_atom_types = self.compute_ligand_atom_counts(ligand)
    hbonds = self.compute_hydrogen_bonds(ligand, receptor)
    hydrophobics = self.compute_hydrophobic_contacts(ligand, receptor)
    stacking = self.compute_pi_pi_stacking(ligand, receptor)
    pi_cation = self.compute_pi_cation(ligand, receptor)
    t_shaped = self.compute_pi_T(ligand, receptor)
    active_site_flexibility = (
      self.compute_active_site_flexibility(ligand, receptor))
    salt_bridges = self.compute_salt_bridges(ligand, receptor)

    input_vector = []
    input_vector = extend_list_by_dictionaries(input_vector, [
        ligand_receptor_contacts,
        ligand_receptor_electrostatics,
        ligand_atom_types,
        ligand_receptor_close_contacts,
        hbonds,
        hydrophobics,
        stacking,
        pi_cation,
        t_shaped,
        active_site_flexibility,
        salt_bridges,
        rotatable_bonds_count])
    return input_vector


  def __init__(self):
    """
    Parameters
    ----------
    """
    self.functions = MathFunctions()
