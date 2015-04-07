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

################################## MODIFY THIS VARIABLE TO POINT TO THE AUTODOCK VINA EXECUTABLE ##################################
# TODO(bramsundar): Hard coding the vina executable is bad design.
# Refactor this to match use of ambertools elsewhere.
#vina_executable = "/PATH/TO/VINA_1_1_2/vina"
vina_executable = (
    "/home/bramsundar/iautodock/autodock_vina_1_1_2_linux_x86/bin/vina")
# Example: vina_executable = "/gpfs/autodock_vina_1_1_2_linux_x86/bin/vina"
###################################################################################################################################


import textwrap
import math
import os
import sys
import textwrap
import glob
import cPickle
from vs_utils.features import Featurizer
from vs_utils.features.nnscore_helper import MathFunctions
from vs_utils.features.nnscore_helper import PDB

# TODO(bramsundar): Many places in this file use long switch
# statements. Could there be a cleaner way to structure this file?

# TODO(bramsundar): Add some tests for the classes and functions in
# this file. Will help verify that cleanup didn't bork anything.

# TODO(bramsundar): Props to the author for the careful featurization,
# but the code structure itself needs a lot of work. Are the misnamed
# variables on purpose in places? I think he's trying to avoid
# variable name collision with misspellings...

# TODO(bramsundar): How does vs_utils handle scripts? Should I factor
# the script part of this file elsewhere?

# TODO(bramsundar): Many arbitrary choices in cutoff-values. I think
# that moving to a 3D-convnet might actually be necessary to remove
# the high-degree of hand-tuning here.

# TODO(bramsundar): Got lazy. Need to reformat some of the help
# strings.

# TODO(bramsundar): The binana feature vectors are organized around
# series of hashes. These hashes aren't described in the comments.
# Fix.

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


class NNScoreFeaturizer(Featurizer):
  """
  Featurizes a given docking pose using the features computed
  by NNScore 2.0.1.
  """

  def _featurize(self, mol):
    """
    Compute NNScore Docking Fingerprint.
    """
    # TODO(bramsundar): Add connector code here.
    pass


#def getCommandOutput2(command):
#  child = os.popen(command)
#  data = child.read()
#  err = child.close()
#  if err:
#      raise RuntimeError, '%s failed w/ exit code %d' % (command, err)
#  return data

def hashtable_entry_add_one(hashtable, key, toadd = 1):
  # note that dictionaries (hashtables) are passed by reference in python
  if hashtable.has_key(key):
    hashtable[key] = hashtable[key] + toadd
  else:
    hashtable[key] = toadd


def extend_list_by_dictionary(lst, dictionary):
  # first, sort the dictionary by the key
  keys = dictionary.keys()
  keys.sort()

  # now make a list of the values
  vals = []
  for key in keys:
    vals.append(dictionary[key])

  # now append vals to the list
  newlist = lst[:]
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

  #functions = MathFunctions()

  # supporting functions

  def get_vina_output(self):
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

  def compute_hydrophobic_contact(self, ligand, receptor):
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
    hydrophobics = {}
    pdb_hydrophobic = PDB()
    for ligand_atom_index in ligand.AllAtoms:
      ligand_atom = ligand.AllAtoms[ligand_atom_index]
      for receptor_atom_index in receptor.AllAtoms:
        receptor_atom = receptor.AllAtoms[receptor_atom_index]
        if ligand_atom.element == "C" and receptor_atom.element == "C":
          hydrophobic_key = (receptor_atom.SideChainOrBackBone() +
              "_" + receptor_atom.structure)
          pdb_hydrophobic.AddNewAtom(ligand_atom.copy_of())
          pdb_hydrophobic.AddNewAtom(receptor_atom.copy_of())

          hashtable_entry_add_one(hydrophobics, hydrophobic_key)
    return hydrophobics, pdb_hydrophobic

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
    ligand_receptor_electrostatics = {}
    for ligand_atom_index in ligand.AllAtoms:
      ligand_atom = ligand.AllAtoms[ligand_atom_index]
      for receptor_atom_index in receptor.AllAtoms:
        receptor_atom = receptor.AllAtoms[receptor_atom_index]
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
    active_site_flexibility = {}
    for receptor_atom_index in receptor.AllAtoms:
      receptor_atom = receptor.AllAtoms[receptor_atom_index]

      flexibility_key = (receptor_atom.SideChainOrBackBone() + "_"
          + receptor_atom.structure)
      hashtable_entry_add_one(active_site_flexibility, flexibility_key)
    return active_site_flexibility


  def compute_hydrogen_bonds(self, ligand, receptor):
    """
    Computes hydrogen bonds between ligand and receptor.

    Returns a dictionary whose keys are of form
    HDONOR_${COMMENT}_${RESIDUETYPE}_${STRUCTURE} where COMMENT is either
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
    hbonds = {}
    for ligand_atom_index in ligand.AllAtoms:
      ligand_atom = ligand.AllAtoms[ligand_atom_index]
      for receptor_atom_index in receptor.AllAtoms:
        receptor_atom = receptor.AllAtoms[receptor_atom_index]
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

          for atm_index in ligand.AllAtoms:
            if ligand.AllAtoms[atm_index].element == "H":
              # so it's a hydrogen
              if (ligand.AllAtoms[atm_index].coordinates.dist_to(
                  ligand_atom.coordinates) < H_BOND_DIST):
                ligand.AllAtoms[atm_index].comment = "LIGAND"
                hydrogens.append(ligand.AllAtoms[atm_index])

          for atm_index in receptor.AllAtoms:
            if receptor.AllAtoms[atm_index].element == "H": # so it's a hydrogen
              if (receptor.AllAtoms[atm_index].coordinates.dist_to(
                  receptor_atom.coordinates) < H_BOND_DIST):
                receptor.AllAtoms[atm_index].comment = "RECEPTOR"
                hydrogens.append(receptor.AllAtoms[atm_index])

          # now we need to check the angles
          # TODO(bramsundar): Verify this heuristic and add a source
          # citation.
          for hydrogen in hydrogens:
            if (math.fabs(180 - functions.angle_between_three_points(
                  ligand_atom.coordinates,
                  hydrogen.coordinates,
                  receptor_atom.coordinates) * 180.0 / math.pi) <=
                  H_BOND_ANGLE):
              hbonds_key = ("HDONOR_" + hydrogen.comment + "_" +
                  receptor_atom.SideChainOrBackBone() + "_" +
                  receptor_atom.structure)
              hashtable_entry_add_one(hbonds, hbonds_key)
    return hbonds

  def compute_ligand_atom_counts(self, ligand):
    """Counts atoms of each type in given ligand.

    Returns a dictionary that makes between atom types ("C", "O", etc.) to
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
    ligand_atom_types = {}
    for ligand_atom_index in ligand.AllAtoms:
      ligand_atom = ligand.AllAtoms[ligand_atom_index]
      hashtable_entry_add_one(ligand_atom_types, ligand_atom.atomtype)
    return ligand_atom_types

  def compute_ligand_receptor_distances_and_charges(self, ligand, receptor):
    """
    Compute distance measurements and electrotstatics between ligand
    and receptor.
    """
    ligand_receptor_close_contacts = {}
    ligand_receptor_contacts = {}
    for ligand_atom_index in ligand.AllAtoms:
      for receptor_atom_index in receptor.AllAtoms:
        ligand_atom = ligand.AllAtoms[ligand_atom_index]
        receptor_atom = receptor.AllAtoms[receptor_atom_index]

        dist = ligand_atom.coordinates.dist_to(receptor_atom.coordinates)
        atomtypes = [ligand_atom.atomtype, receptor_atom.atomtype]
        atomstr = "_".join(sorted(atomtypes))
        if dist < CLOSE_CONTACT_CUTOFF:
          hashtable_entry_add_one(ligand_receptor_close_contacts,
              atomstr)
          pdb_close_contacts.AddNewAtom(ligand_atom.copy_of())
          pdb_close_contacts.AddNewAtom(receptor_atom.copy_of())
        elif dist < CONTACT_CUTOFF:
          hashtable_entry_add_one(ligand_receptor_contacts,
              atomstr)

          #self.compute_electrostatic_energy(ligand_atom, receptor_atom)
          #self.compute_active_site_flexibility(ligand_atom, receptor_atom)
          #self.compute_hydrophobic_contact(ligand_atom, receptor_atom)
          #self.compute_hydrogen_bond(ligand_atom, receptor_atom)

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
    for lig_aromatic in ligand.aromatic_rings:
      for rec_aromatic in receptor.aromatic_rings:
        dist = lig_aromatic.center.dist_to(rec_aromatic.center)
        if dist < PI_PI_CUTOFF:
          # so there could be some pi-pi interactions.  Now, let's
          # check for stacking interactions. Are the two pi's roughly
          # parallel?
          lig_aromatic_norm_vector = point(lig_aromatic.plane_coeff[0],
              lig_aromatic.plane_coeff[1], lig_aromatic.plane_coeff[2])
          rec_aromatic_norm_vector = point(rec_aromatic.plane_coeff[0],
              rec_aromatic.plane_coeff[1], rec_aromatic.plane_coeff[2])
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
                  ligand.AllAtoms[ligand_ring_index].coordinates,
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
                    receptor.AllAtoms[receptor_ring_index].coordinates,
                    lig_aromatic.plane_coeff)
                if (pt_on_ligand_plane.dist_to(lig_aromatic.center)
                    <= lig_aromatic.radius + PI_PADDING):
                  pi_pi = True
                  break

            if pi_pi == True:
                structure = receptor.AllAtoms[rec_aromatic.indices[0]].structure
                if structure == "":
                  # since it could be interacting with a cofactor or something
                  structure = "OTHER"
                key = "STACKING_" + structure


                hashtable_entry_add_one(PI_interactions, key)

          elif math.fabs(angle_between_planes-90) < 30.0 or math.fabs(angle_between_planes-270) < 30.0:
            # so they're more or less perpendicular, it's probably a
            # pi-edge interaction having looked at many structures, I
            # noticed the algorithm was identifying T-pi reactions
            # when the two rings were in fact quite distant, often
            # with other atoms in between. Eye-balling it, requiring
            # that at their closest they be at least 5 A apart seems
            # to separate the good T's from the bad
            min_dist = 100.0
            for ligand_ind in lig_aromatic.indices:
              ligand_at = ligand.AllAtoms[ligand_ind]
              for receptor_ind in rec_aromatic.indices:
                receptor_at = receptor.AllAtoms[receptor_ind]
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
                structure = receptor.AllAtoms[rec_aromatic.indices[0]].structure
                if structure == "":
                  # since it could be interacting with a cofactor or something
                  structure = "OTHER"
                key = "T-SHAPED_" + structure

                hashtable_entry_add_one(PI_interactions, key)

  def compute_pi_cation(self, ligand, receptor):
    """
    Computes number of pi-cation interactions.

    Parameters
    ----------
    ligand: PDB Object
      small molecule to dock.
    receptor: PDB Object
      protein to dock agains.
    """
    for aromatic in receptor.aromatic_rings:
      for charged in ligand.charges:
        if charged.positive == True: # so only consider positive charges
          if (charged.coordinates.dist_to(aromatic.center) < CATION_PI_CUTOFF):

            # project the charged onto the plane of the aromatic
            charge_projected = self.functions.project_point_onto_plane(
                charged.coordinates,aromatic.plane_coeff)
            if (charge_projected.dist_to(aromatic.center)
              < aromatic.radius + PI_PADDING):
              structure = receptor.AllAtoms[aromatic.indices[0]].structure
              if structure == "":
                # since it could be interacting with a cofactor or something
                structure = "OTHER"
              key = "PI-CATION_LIGAND-CHARGED_" + structure

              for index in aromatic.indices:
                pdb_pi_cat.AddNewAtom(receptor.AllAtoms[index].copy_of())
              for index in charged.indices:
                pdb_pi_cat.AddNewAtom(ligand.AllAtoms[index].copy_of())

              hashtable_entry_add_one(PI_interactions, key)

    for aromatic in ligand.aromatic_rings:
      # now it's the ligand that has the aromatic group
      for charged in receptor.charges:
        if charged.positive == True: # so only consider positive charges
          if (charged.coordinates.dist_to(aromatic.center) < CATION_PI_CUTOFF):
            charge_projected = self.functions.project_point_onto_plane(
                charged.coordinates,aromatic.plane_coeff)
            if charge_projected.dist_to(aromatic.center) < aromatic.radius + PI_PADDING:
              structure = receptor.AllAtoms[charged.indices[0]].structure
              if structure == "":
                # since it could be interacting with a cofactor or something
                structure = "OTHER"
              key = "PI-CATION_RECEPTOR-CHARGED_" + structure

              for index in aromatic.indices:
                pdb_pi_cat.AddNewAtom(ligand.AllAtoms[index].copy_of())
              for index in charged.indices:
                pdb_pi_cat.AddNewAtom(receptor.AllAtoms[index].copy_of())

              hashtable_entry_add_one(PI_interactions, key)

  def compute_salt_bridges(self, ligand, receptor):
    """
    Computes number of ligand-receptor salt bridges.

    Parameters
    ----------
    ligand: PDB Object
      small molecule to dock.
    receptor: PDB Object
      protein to dock agains.
    """
    for receptor_charge in receptor.charges:
      for ligand_charge in ligand.charges:
        if ligand_charge.positive != receptor_charge.positive:
          # so they have oppositve charges
          if (ligand_charge.coordinates.dist_to(
              receptor_charge.coordinates) < SALT_BRIDGE_CUTOFF):
            structure = receptor.AllAtoms[receptor_charge.indices[0]].structure
            if structure == "":
              # since it could be interacting with a cofactor or something
              structure = "OTHER"
            key = "SALT-BRIDGE_" + structure

            for index in receptor_charge.indices:
              pdb_salt_bridges.AddNewAtom(receptor.AllAtoms[index].copy_of())
            for index in ligand_charge.indices:
              pdb_salt_bridges.AddNewAtom(ligand.AllAtoms[index].copy_of())

            hashtable_entry_add_one(salt_bridges, key)

  def compute_data(self):
    # now create a single descriptor object
    data = {}
    data['vina_output'] = self.get_vina_output()
    data['ligand_receptor_close_contacts'] = (
        ligand_receptor_close_contacts)
    data['ligand_receptor_contacts'] = (
        ligand_receptor_contacts)
    data['ligand_atom_types'] = ligand_atom_types
    data['ligand_receptor_electrostatics'] = (
        ligand_receptor_electrostatics)
    data['rotatable_bonds_count'] = ligand.rotatable_bonds_count
    data['active_site_flexibility'] = active_site_flexibility
    data['hbonds'] = hbonds
    data['hydrophobics'] = hydrophobics
    data['stacking'] = stacking
    data['pi_cation'] = pi_cation
    data['t_shaped'] = t_shaped
    data['salt_bridges'] = salt_bridges
    return data

  def compute_input_vector(self):
    # TODO(bramsundar): There is a ton of data copying in this
    # function? Can I strip this out? It looks like a bunch of these
    # terms are not used anywhere else in the file. Might be ok to
    # just delete it. Ok, I do see one use. It looks like this step is
    # used to provide padding so we have a fixed-length vector.
    data = self.compute_data()
    vina_output = data['vina_output']
    rotatable_bonds_count = {'rot_bonds': data['rotatable_bonds_count']}

    active_site_flexibility = {
      'BACKBONE_ALPHA': 0, 'BACKBONE_BETA': 0, 'BACKBONE_OTHER': 0,
      'SIDECHAIN_ALPHA': 0, 'SIDECHAIN_BETA': 0, 'SIDECHAIN_OTHER': 0
      }
    for key in data['active_site_flexibility']:
      active_site_flexibility[key] = data['active_site_flexibility'][key]

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
    for key in data['hbonds']:
      key2 = key.replace("HDONOR_","HDONOR-")
      hbonds[key2] = data['hbonds'][key]


    hydrophobics = {
      'BACKBONE_ALPHA': 0, 'BACKBONE_BETA': 0, 'BACKBONE_OTHER': 0,
      'SIDECHAIN_ALPHA': 0, 'SIDECHAIN_BETA': 0, 'SIDECHAIN_OTHER': 0
      }
    for key in data['hydrophobics']:
      hydrophobics[key] = data['hydrophobics'][key]


    stacking_tmp = {}
    for item in data['stacking']:
      item = item.split("_")
      stacking_tmp[item[1]] = int(item[2])
    stacking = {'ALPHA': 0, 'BETA': 0, 'OTHER': 0}
    for key in stacking_tmp:
      stacking[key] = stacking_tmp[key]

    all = stacking['ALPHA'] + stacking['BETA'] + stacking['OTHER']
    stacking_all = {'all': all}

    pi_cation_tmp = {}
    for item in data['pi_cation']:
      item = item.split("_")
      pi_cation_tmp[item[1] + "_" + item[2]] = int(item[3])
    pi_cation = {
      'LIGAND-CHARGED_ALPHA': 0, 'LIGAND-CHARGED_BETA': 0,
      'LIGAND-CHARGED_OTHER': 0, 'RECEPTOR-CHARGED_ALPHA': 0,
      'RECEPTOR-CHARGED_BETA': 0, 'RECEPTOR-CHARGED_OTHER': 0}
    for key in pi_cation_tmp:
      pi_cation[key] = pi_cation_tmp[key]


    t_shaped_tmp = {}
    for item in data['t_shaped']:
      item = item.split("_")
      t_shaped_tmp[item[1]] = int(item[2])
    t_shaped = {'ALPHA': 0, 'BETA': 0, 'OTHER': 0}
    for key in t_shaped_tmp:
      t_shaped[key] = t_shaped_tmp[key]

    salt_bridges = {'ALPHA': 0, 'BETA': 0, 'OTHER': 0}
    for key in data['salt_bridges']:
      key2 = key.replace("SALT-BRIDGE_","")
      salt_bridges[key2] = data['salt_bridges'][key]

    input_vector = []
    input_vector.extend(vina_output) # a list
    input_vector = extend_list_by_dictionary(
        input_vector, ligand_receptor_contacts)
    input_vector = extend_list_by_dictionary(
        input_vector, ligand_receptor_electrostatics)
    input_vector = extend_list_by_dictionary(
        input_vector, ligand_atom_types)
    input_vector = extend_list_by_dictionary(
        input_vector, ligand_receptor_close_contacts)
    input_vector = extend_list_by_dictionary(
        input_vector, hbonds)
    input_vector = extend_list_by_dictionary(
        input_vector, hydrophobics)
    input_vector = extend_list_by_dictionary(
        input_vector, stacking)
    input_vector = extend_list_by_dictionary(
        input_vector, pi_cation)
    input_vector = extend_list_by_dictionary(
        input_vector, t_shaped)
    input_vector = extend_list_by_dictionary(
        input_vector, active_site_flexibility)
    input_vector = extend_list_by_dictionary(
        input_vector, salt_bridges)
    input_vector = extend_list_by_dictionary(
        input_vector, rotatable_bonds_count)
    return input_vector


  def process_pi(self):
    stacking = []
    t_shaped = []
    pi_cation = []
    for key in PI_interactions:
      value = PI_interactions[key]
      together = key + "_" + str(value) # not that this is put together strangely!!!
      if "STACKING" in together: stacking.append(together)
      if "CATION" in together: pi_cation.append(together)
      if "SHAPED" in together: t_shaped.append(together)


  # The meat of the class; must be a more elegant way of doing this
  # TODO(bramsundar): Factor this __init__ into multiple functions to
  # avoid ~ 800 line function.
  #def __init__(self, ligand_pdbqt_filename, receptor, parameters, line_header):
  def __init__(self):
    """
    Parameters
    ----------
    ligand_pdbqt_filename: string
      path to ligand's pdbqt file (generated by vina)
    receptor: Object of type PDB for the receptor
    parameters: Object of type command_line_parameters.
      TODO(bramsundar): Having a featurizer accept a
      command_line_parameters object seems like a failure of
      abstraction.
    line_header: string
      TODO(bramsundar): line separator in PDB files (?)
    """

    pass
#   receptor_pdbqt_filename = receptor.OrigFileName
#
#   self.ligand = PDB()
#   self.ligand.LoadPDB_from_file(ligand_pdbqt_filename, line_header)
#
#   receptor.assign_secondary_structure()

#   # Get distance measurements between protein and ligand atom types,
#   # as well as some other measurements
#   ligand.rotatable_bonds_count
#   functions = MathFunctions()

#   # These dictionaries are keyed by strings of atom-pairs (C_O for
#   # example)
#   ligand_receptor_close_contacts = {}
#   ligand_receptor_contacts = {}
#
#   pdb_close_contacts = PDB()

#   self.compute_ligand_receptor_distances_and_charges()
#
#   # Get the total number of each atom type in the ligand
#   ligand_atom_types = self.compute_ligand_atom_counts(
#       self.ligand)
#
#
#
#   # Count pi-pi stacking and pi-T stacking interactions
#   self.PI_interactions = {}

#   self.compute_pi_pi_stacking(self.ligand, self.receptor)
#
#   # Now identify pi-cation interactions
#   self.compute_pi_cation(self.ligand, self.receptor)
#   self.pdb_pi_cat = PDB()
#
#
#   # now count the number of salt bridges
#   self.pdb_salt_bridges = PDB()
#   self.salt_bridges = {}

#
#
#
#    # TODO(bramsundar): Move this elsewhere?
#    self.ligand_receptor_close_contacts = {
#        "A_A": 0, "A_C": 0, "A_CL": 0, "A_F": 0, "A_FE": 0, "A_MG": 0,
#        "A_MN": 0, "A_NA": 0, "A_SA": 0, "BR_C": 0, "BR_OA": 0, "C_CL":
#        0, "CD_OA": 0, "CL_FE": 0, "CL_MG": 0, "CL_N": 0, "CL_OA": 0,
#        "CL_ZN": 0, "C_MN": 0, "C_NA": 0, "F_N": 0, "F_SA": 0, "F_ZN":
#        0, "HD_MN": 0, "MN_N": 0, "NA_SA": 0, "N_SA": 0, "A_HD": 0,
#        "A_N": 0, "A_OA": 0, "A_ZN": 0, "BR_HD": 0, "C_C": 0, "C_F": 0,
#        "C_HD": 0, "CL_HD": 0, "C_MG": 0, "C_N": 0, "C_OA": 0, "C_SA":
#        0, "C_ZN": 0, "FE_HD": 0, "FE_N": 0, "FE_OA": 0, "F_HD": 0,
#        "F_OA": 0, "HD_HD": 0, "HD_I": 0, "HD_MG": 0, "HD_N": 0,
#        "HD_NA": 0, "HD_OA": 0, "HD_P": 0, "HD_S": 0, "HD_SA": 0,
#        "HD_ZN": 0, "MG_NA": 0, "MG_OA": 0, "MN_OA": 0, "NA_OA": 0,
#        "NA_ZN": 0, "N_N": 0, "N_NA": 0, "N_OA": 0, "N_ZN": 0, "OA_OA":
#        0, "OA_SA": 0, "OA_ZN": 0, "SA_ZN": 0, "S_ZN": 0}
#    for key in data['ligand_receptor_close_contacts']:
#      if not key in self.ligand_receptor_close_contacts:
#          print "\tWARNING: Atoms of types " + key.replace("_"," and ") + " come within 2.5 angstroms of each other."
#          print "\t  The neural networks were not trained to deal with this juxtaposition,"
#          print "\t  so it will be ignored."
#          self.error = True
#      else:
#         self.ligand_receptor_close_contacts[key] = (
#             data['ligand_receptor_close_contacts'][key])
#
#    self.ligand_receptor_contacts = {
#        "A_CU": 0, "A_MG": 0, "A_MN": 0, "BR_SA": 0, "C_CD": 0,
#        "CL_FE": 0, "CL_MG": 0, "CL_MN": 0, "CL_NA": 0, "CL_P": 0,
#        "CL_S": 0, "CL_ZN": 0, "CU_HD": 0, "CU_N": 0, "FE_NA": 0,
#        "FE_SA": 0, "MG_N": 0, "MG_S": 0, "MG_SA": 0, "MN_NA": 0,
#        "MN_S": 0, "MN_SA": 0, "NA_P": 0, "P_S": 0, "P_SA": 0, "S_SA":
#        0, "A_A": 0, "A_BR": 0, "A_C": 0, "A_CL": 0, "A_F": 0, "A_FE":
#        0, "A_HD": 0, "A_I": 0, "A_N": 0, "A_NA": 0, "A_OA": 0, "A_P":
#        0, "A_S": 0, "A_SA": 0, "A_ZN": 0, "BR_C": 0, "BR_HD": 0,
#        "BR_N": 0, "BR_OA": 0, "C_C": 0, "C_CL": 0, "C_F": 0, "C_FE":
#        0, "C_HD": 0, "C_I": 0, "CL_HD": 0, "CL_N": 0, "CL_OA": 0,
#        "CL_SA": 0, "C_MG": 0, "C_MN": 0, "C_N": 0, "C_NA": 0, "C_OA":
#        0, "C_P": 0, "C_S": 0, "C_SA": 0, "C_ZN": 0, "FE_HD": 0,
#        "FE_N": 0, "FE_OA": 0, "F_HD": 0, "F_N": 0, "F_OA": 0, "F_SA":
#        0, "HD_HD": 0, "HD_I": 0, "HD_MG": 0, "HD_MN": 0, "HD_N": 0,
#        "HD_NA": 0, "HD_OA": 0, "HD_P": 0, "HD_S": 0, "HD_SA": 0,
#        "HD_ZN": 0, "I_N": 0, "I_OA": 0, "MG_NA": 0, "MG_OA": 0,
#        "MG_P": 0, "MN_N": 0, "MN_OA": 0, "MN_P": 0, "NA_OA": 0,
#        "NA_S": 0, "NA_SA": 0, "NA_ZN": 0, "N_N": 0, "N_NA": 0,
#        "N_OA": 0, "N_P": 0, "N_S": 0, "N_SA": 0, "N_ZN": 0, "OA_OA":
#        0, "OA_P": 0, "OA_S": 0, "OA_SA": 0, "OA_ZN": 0, "P_ZN": 0,
#        "SA_SA": 0, "SA_ZN": 0, "S_ZN": 0}
#    for key in data['ligand_receptor_contacts']:
#      if not key in self.ligand_receptor_contacts:
#        print "\tWARNING: Atoms of types " + key.replace("_"," and ") + " come within 4 angstroms of each other."
#        print "\t  The neural networks were not trained to deal with this juxtaposition,"
#        print "\t  so it will be ignored."
#
#        self.error = True
#      else:
#        self.ligand_receptor_contacts[key] = (
#            data['ligand_receptor_contacts'][key])
#
#    self.ligand_atom_types = {
#        'A': 0, 'BR': 0, 'C': 0, 'CL': 0, 'F': 0, 'HD': 0, 'I': 0,
#        'N': 0, 'NA': 0, 'OA': 0, 'P': 0, 'S': 0, 'SA': 0}
#    for key in data['ligand_atom_types']:
#      if not key in self.ligand_atom_types:
#          print "\tWARNING: The ligand contains an atoms of type " + key + ". The neural networks"
#          print "\t  were not trained to deal with this ligand atom type, so it will be ignored."
#
#          self.error = True
#      else:
#          self.ligand_atom_types[key] = data['ligand_atom_types'][key]
#
#    self.ligand_receptor_electrostatics = {
#        "A_MG": 0, "A_MN": 0, "BR_SA": 0, "CL_FE": 0, "CL_MG": 0,
#        "CL_MN": 0, "CL_NA": 0, "CL_P": 0, "CL_S": 0, "CL_ZN": 0,
#        "CU_HD": 0, "CU_N": 0, "FE_NA": 0, "FE_SA": 0, "MG_N": 0,
#        "MG_S": 0, "MG_SA": 0, "MN_NA": 0, "MN_S": 0, "MN_SA": 0,
#        "NA_P": 0, "P_S": 0, "P_SA": 0, "S_SA": 0, "A_A": 0.0, "A_BR":
#        0.0, "A_C": 0.0, "A_CL": 0.0, "A_F": 0.0, "A_FE": 0.0, "A_HD":
#        0.0, "A_I": 0.0, "A_N": 0.0, "A_NA": 0.0, "A_OA": 0.0, "A_P":
#        0.0, "A_S": 0.0, "A_SA": 0.0, "A_ZN": 0.0, "BR_C": 0.0,
#        "BR_HD": 0.0, "BR_N": 0.0, "BR_OA": 0.0, "C_C": 0.0, "C_CL":
#        0.0, "C_F": 0.0, "C_FE": 0.0, "C_HD": 0.0, "C_I": 0.0,
#        "CL_HD": 0.0, "CL_N": 0.0, "CL_OA": 0.0, "CL_SA": 0.0, "C_MG":
#        0.0, "C_MN": 0.0, "C_N": 0.0, "C_NA": 0.0, "C_OA": 0.0, "C_P":
#        0.0, "C_S": 0.0, "C_SA": 0.0, "C_ZN": 0.0, "FE_HD": 0.0,
#        "FE_N": 0.0, "FE_OA": 0.0, "F_HD": 0.0, "F_N": 0.0, "F_OA":
#        0.0, "F_SA": 0.0, "HD_HD": 0.0, "HD_I": 0.0, "HD_MG": 0.0,
#        "HD_MN": 0.0, "HD_N": 0.0, "HD_NA": 0.0, "HD_OA": 0.0, "HD_P":
#        0.0, "HD_S": 0.0, "HD_SA": 0.0, "HD_ZN": 0.0, "I_N": 0.0,
#        "I_OA": 0.0, "MG_NA": 0.0, "MG_OA": 0.0, "MG_P": 0.0, "MN_N":
#        0.0, "MN_OA": 0.0, "MN_P": 0.0, "NA_OA": 0.0, "NA_S": 0.0,
#        "NA_SA": 0.0, "NA_ZN": 0.0, "N_N": 0.0, "N_NA": 0.0, "N_OA":
#        0.0, "N_P": 0.0, "N_S": 0.0, "N_SA": 0.0, "N_ZN": 0.0,
#        "OA_OA": 0.0, "OA_P": 0.0, "OA_S": 0.0, "OA_SA": 0.0, "OA_ZN":
#        0.0, "P_ZN": 0.0, "SA_SA": 0.0, "SA_ZN": 0.0, "S_ZN": 0,
#        "F_ZN": 0}
#
#    for key in data['ligand_receptor_electrostatics']:
#      if not key in self.ligand_receptor_electrostatics:
#         print "\tWARNING: Atoms of types " + key.replace("_"," and ") + ", which come within 4 angstroms of each"
#         print "\t  other, may interact electrostatically. However, the neural networks"
#         print "\t  were not trained to deal with electrostatic interactions between atoms"
#         print "\t  of these types, so they will be ignored."
#
#         self.error = True
#      else:
#         self.ligand_receptor_electrostatics[key] = (
#             data['ligand_receptor_electrostatics'][key])
#


# TODO(bramsundar): I think vs_utils has a standard way of dealing
# with command_line parameters. Check scripts.
class command_line_parameters:

  params = {}

  def __init__(self, parameters):

    global vina_executable
#
#    # first, set defaults
#    self.params['receptor'] = ''
#    self.params['ligand'] = ''
#    self.params['vina_executable'] = vina_executable
#    # TRUE by default, but setting to false will speed up execution.
#    # Good when rescoring many poses.
#    self.params['check_vina_version'] = "TRUE"
#
#    # now get user inputed values
#
#    for index in range(len(parameters)):
#      item = parameters[index]
#      if item[:1] == '-': # so it's a parameter key value
#        key = item.replace('-','').lower()
#
#        if key == "help":
#          print "INTRODUCTION"
#          print "============\n"
#          print textwrap.fill("NNScore 2.01 (NNScore2.py) is a "
#          "python script for predicting the binding affinity of "
#          "receptor-ligand complexes. It is particularly well "
#          "suited for rescoring ligands that have been docked
#          using AutoDock Vina.") + "\n"
#          print "REQUIREMENTS"
#          print "============\n"
#          print textwrap.fill("Python: NNScore 2.01 has been "
#          "tested using Python versions 2.6.5, 2.6.1, and 2.5.2 on "
#          "Ubuntu 10.04.1 LTS, Mac OS X 10.6.8, and Windows XP "
#          "Professional, respectively. A copy of the Python "
#          "interpreter can be downloaded from "
#          "http://www.python.org/getit/.") + "\n"
#          print textwrap.fill("AutoDock Vina 1.1.2: NNScore 2.01 "
#          "uses AutoDock Vina 1.1.2 to obtain some information "
#          "about the receptor-ligand complex. Note that previous "
#          "versions of AutoDock Vina are not suitble. AutoDock Vina "
#          "1.1.2 can be downloaded from "
#          "http://vina.scripps.edu/download.html.") + "\n"
#          print textwrap.fill("MGLTools: As receptor and ligand
#          inputs, NNScore 2.01 accepts models in the PDBQT format.
#          Files in the more common PDB format can be converted to
#          the PDBQT format using scripts included in MGLTools
#          (prepare_receptor4.py and prepare_ligand4.py). MGLTools
#          can be obtained from
#          http://mgltools.scripps.edu/downloads.") + "\n"
#          print "COMMAND-LINE PARAMETERS"
#          print "=======================\n"
#          print textwrap.fill("-receptor: File name of the
#          receptor PDBQT file.") + "\n"
#          print textwrap.fill("-ligand: File name of the ligand
#          PDBQT file. AutoDock Vina output files, typically
#          containing multiple poses, are also permitted.") + "\n"
#          print textwrap.fill("-vina_executable: The location of
#          the AutoDock Vina 1.1.2 executable. If you don't wish to
#          specify the location of this file every time you use
#          NNScore 2.01, simply edit the vina_executable variable
#          defined near the begining of the NNScore2.py script.") +
#          "\n"
#          print "PROGRAM OUTPUT"
#          print "==============\n"
#          print textwrap.fill("NNScore 2.01 evaluates each of the
#          ligand poses contained in the file specified by the
#          -ligand tag using 20 distinct neural-network scoring
#          functions. The program then seeks to identify which of
#          the poses has the highest predicted affinity using
#          several metrics:") + "\n"
#          print textwrap.fill("1) Each of the 20 networks are
#          considered separately. The poses are ranked in 20
#          different ways by the scores assigned by each network.") + "\n"
#          print textwrap.fill("2) The poses are ranked by the best
#          score given by any of the 20 networks.") + "\n"
#          print textwrap.fill("3) The poses are ranked by the
#          average of the scores given by the 20 networks. This is
#          the recommended metric.") + "\n"
#
#
#          sys.exit(0)
#
#        value = parameters[index+1]
#        if key in self.params.keys():
#            self.params[key] = value
#            parameters[index] = ""
#            parameters[index + 1] = ""
#
#    if self.okay_to_proceed() == True:
#      print "Command-line parameters used:"
#      print "\tReceptor:        " + self.params["receptor"]
#      print "\tLigand:          " + self.params["ligand"]
#      print "\tVina executable: " + self.params["vina_executable"] + "\n"
#
#    # make a list of all the command-line parameters not used
#    error = ""
#    for index in range(1,len(parameters)):
#      item = parameters[index]
#      if item != "": error = error + item + " "
#
#    if error != "":
#      print "WARNING: The following command-line parameters were not used:"
#      print "\t" + error + "\n"
#
#
#    # do a check to make sure the autodock vina is version 1.1.2
#    if self.params["check_vina_version"] == "TRUE":
#      vina_version_output = ""
#      if not os.path.exists(self.params['vina_executable']):
#        vina_version_output = ""
#      else:
#        try:
#          vina_version_output = getCommandOutput2(self.params['vina_executable'] + ' --version')
#        except:
#          vina_version_output = ""
#
#      if not " 1.1.2 " in vina_version_output:
#        print "ERROR: NNScore 2.01 is designed to work with AutoDock
#        Vina 1.1.2.\nOther versions of AutoDock may not work properly.
#        Please download\nAutoDock Vina 1.1.2 from
#        http://vina.scripps.edu/download.html.\n"
#        print "Once this executable is downloaded, you can use the
#        -vina_executable\ntag to indicate its location. Alternatively,
#        you can simply modify\nthe vina_executable variable defined
#        near the beginning of\nthe NNScore2.py file.\n"
#        sys.exit(0)
#
#
#  def okay_to_proceed(self):
#      if self.params['receptor'] != '' and self.params['ligand'] != '' and self.params['vina_executable'] != '':
#          return True
#      else: return False

#def score_to_kd(score):
#  kd = math.pow(10,-score)
#  if score <= 0: return "Kd = " + str(round(kd,2)) + " M"
#  temp = kd * pow(10,3)
#  if temp >=1 and temp <1000: return "Kd = " + str(round(temp,2)) + " mM"
#  temp = temp * pow(10,3)
#  if temp >=1 and temp <1000: return "Kd = " + str(round(temp,2)) + " uM"
#  temp = temp * pow(10,3)
#  if temp >=1 and temp <1000: return "Kd = " + str(round(temp,2)) + " nM"
#  temp = temp * pow(10,3)
#  if temp >=1 and temp <1000: return "Kd = " + str(round(temp,2)) + " pM"
#  temp = temp * pow(10,3)
#  #if temp >=1 and temp <1000:
#  return "Kd = " + str(round(temp,2)) + " fM"
#  #return "Kd = " + str(kd) + " M"

#print "\nNNScore 2.01"
#print "============\n"
#
#print textwrap.fill("NNScore 2.01 is released under the GNU General
#Public License (see http://www.gnu.org/licenses/gpl.html). If you have
#any questions, comments, or suggestions, please contact the author,
#Jacob Durrant, at jdurrant [at] ucsd [dot] edu. If you use NNScore
#2.01 in your work, please cite [REFERENCE HERE].")
#
#print ""
#
#print textwrap.fill("NNScore 2.01 is based in part on the ffnet python
#package developed by Marek Wojciechowski
#(http://ffnet.sourceforge.net/).")
#print ""
#
#print "Use the -help command-line parameter for extended help.\n"
#print "Example: python NNScore2.py -receptor myreceptor.pdbqt -ligand
#myligand.pdbqt -vina_executable /PATH/TO/VINA/1.1.2/vina\n"
#
#cmd_params = command_line_parameters(sys.argv[:])
#
#if cmd_params.okay_to_proceed() == False:
#    print "ERROR: You need to specify the ligand and receptor PDBQT
#    files, as\nwell as the full path the the AutoDock Vina 1.1.2
#    executable, using the\n-receptor, -ligand, and -vina_executable
#    tags from the command line.\nThe -ligand tag can also specify an
#    AutoDock Vina output file.\n"
#    sys.exit(0)
#
#lig = cmd_params.params['ligand']
#rec = cmd_params.params['receptor']

#def calculate_score(lig, rec, cmd_params,
#    actual_filename_if_lig_is_list="", actual_filename_if_rec_is_list="",
#    line_header = ""):
#
#  d = binana(lig, rec, cmd_params, line_header,
#      actual_filename_if_lig_is_list, actual_filename_if_rec_is_list)
#
#  # now load in neural networks
#  scores = []
#  total = 0.0
#  nets = networks()
#  for net_array in nets:
#    try:
#      net = ffnet()
#      net.load(net_array)
#
#      val = net.normcall(d.input_vector)
#      print (line_header + "Network #" + str(len(scores) + 1) + " gave
#          a score of " + str( round(val,3)) + " (" + score_to_kd(val) +
#          ")")
#      scores.append(val)
#      total = total + val
#    except OverflowError:
#      print (line_header + "The output of network #" +
#          str(len(scores)) + 1 + " could not be determined because of an
#          overflow error!" )
#
#  if len(scores) != 0:
#    average = total / len(scores)
#
#     best_score = 0.0
#     best_network = -1
#     count = 1
#    sum = 0.0
#    for score in scores:
#      if score > best_score:
#          best_score = score
#          best_network = count
#          sum = sum + math.pow(score - average,2)
#                      count = count + 1
#    stdev = math.pow(sum / (len(scores)-1), 0.5)
#
#    print ""
#
#    print line_header + "Best Score:         ", round(best_score,3), "(" + score_to_kd(best_score) + ")"
#    print line_header + "Average Score:      ", round(average,3), "(" + score_to_kd(average) + ")"
#    print line_header + "Standard Deviation: ", round(stdev,3)
#    print ""
#
#    return [average, stdev, best_score, best_network, scores]
#  else:
#    print (line_header + "Could not compute the score of this
#    receptor-ligand complex because none of the networks returned a
#    valid score.")
#    return [0, 0, 0, 0, []]

# load the rec into an array so you only have to load it from the disk once
#print "\nLOADING THE RECEPTOR"
#print "====================\n"
#
#receptor = PDB()
#receptor.LoadPDB_from_file(rec)
#receptor.OrigFileName = rec
#
#print "\nEVALUATING EACH OF THE POSES IN THE LIGAND FILE USING 20 TRAINED NEURAL NETWORKS"
#print "================================================================================\n"
#
## determine if the ligand input file is a single pdbqt or an autodock
## vina output file. Both are acceptable inputs.
#
#f = open(lig, 'r')
#
#lig_array = []
#line = "NULL"
#
#scores = []
#
#model_id = 1
#while len(line) != 0:
#  line = f.readline()
#  if line[:6] != "ENDMDL": lig_array.append(line)
#  if line[:6] == "ENDMDL" or len(line) == 0:
#    if len(lig_array) != 0 and lig_array != ['']:
#      temp_filename = lig + ".MODEL_" + str(model_id) + ".pdbqt"
#
#      temp_f = open(temp_filename, 'w')
#      for ln in lig_array: temp_f.write(ln)
#      temp_f.close()
#
#      model_name = "MODEL " + str(model_id)
#
#      print model_name
#
#      score=calculate_score(lig_array, receptor, cmd_params, temp_filename, rec, "\t")
#      scores.append([score[0], score[1], score[2], score[3], score[4], model_name])
#
#      os.remove(temp_filename)
#
#      lig_array = []
#      model_id = model_id + 1
#
#f.close()
#
## now find the best score across all models, for each network
#print "\nRANKED POSES AND SCORES WHEN EACH OF THE 20 NETWORKS IS CONSIDERED SEPARATELY"
#print "=============================================================================\n"
#
#best_network_scores = []
#for t in range(20):
#  net_scores = []
#  for score in scores:
#      net_scores.append((score[4][t],score[5]))
#  net_scores = sorted(net_scores,key=lambda net_score: net_score[0], reverse=True) # sort by score
#  print ("USING NETWORK #" + str(t+1)).center(42," ")
#  print " Rank | Pose     | Score | Predicted Kd "
#  print "------+----------+-------+--------------"
#  count = 1
#  best_network_scores.append(net_scores[0])
#  for net_score in net_scores:
#    print str(count).center(6," ") + "|" + net_score[1].center(10," ") + "|" + str(round(net_score[0],3)).center(7," ") + "|" + score_to_kd(net_score[0]).replace("Kd = ","").center(14," ")
#      count = count + 1
#  print ""
#
#print "\nRANKED POSES AND SCORES WHEN A CONSENSUS OF NETWORK OUTPUTS ARE CONSIDERED"
#print "==========================================================================\n"
#
## now get the best score of all vina output poses
#best_scores = sorted(scores,key=lambda score: score[2], reverse=True) # sort by score
#
#best_best_score = best_scores[0]
#
#print "BEST SCORE OF ALL 20 NETWORKS, BY POSE".center(65," ")
#print " Rank | Pose     | Average Score | Predicted Kd | Network"
#print "------+----------+---------------+--------------+---------"
#
#count = 1
#for score in best_scores:
#  print str(count).center(6," ") + "|" + score[5].center(10," ") + "|" + str(round(score[2],3)).center(15," ") + "|" + score_to_kd(score[2]).replace("Kd = ","").center(14," ") + "|" + ("#" + str(score[3])).center(9," ")
#  count = count + 1
#
#print ""
#
## now get the best average score of all vina output poses
#average_scores = sorted(scores,key=lambda score: score[0], reverse=True) # sort by score
#
#best_of_average_scores = average_scores[0]
#
#print "AVERAGE SCORE OF ALL 20 NETWORKS, BY POSE".center(70," ")
#print " Rank | Pose     | Average Score | Standard Deviation | Predicted Kd"
#print "------+----------+---------------+--------------------+--------------"
#
#count = 1
#for score in average_scores:
#  print str(count).center(6," ") + "|" + score[5].center(10," ") + "|" + str(round(score[0],3)).center(15," ") + "|" + str(round(score[1],3)).center(20," ") + "|" + score_to_kd(score[0]).replace("Kd = ","").center(15," ")
#  count = count + 1
#
#print ""
#print "\nSUMMARY"
#print "=======\n"
#
#count = 1
#for best_network_score in best_network_scores:
#  print textwrap.fill("Best pose scored by network #" + str(count) + ": " + best_network_score[1] + " (Score = " + str(round(best_network_score[0],3)) + " = " + score_to_kd(best_network_score[0]).replace("Kd = ","") + ")")
#  count = count + 1
#
#print ""
#
#print textwrap.fill("When the poses were ranked by the best of the 20 network scores associated with each pose, the best-scoring pose was " + best_best_score[5] + " (Score = " + str(round(best_best_score[2],3)) + " = " + score_to_kd(best_best_score[2]).replace("Kd = ","") + ")")
#print ""
#print textwrap.fill("When the poses were ranked by the average of the 20 network scores associated with each pose, the best-scoring pose was " + best_of_average_scores[5] + " (Score = " + str(round(best_of_average_scores[0],3)) + " +/- " + str(round(best_of_average_scores[1],3)) + " = " + score_to_kd(best_of_average_scores[0]).replace("Kd = ","") + "). This is the recommended ranking/scoring metric.")
#
#print ""
