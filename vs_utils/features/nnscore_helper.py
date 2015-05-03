"""
Helper Classes and Functions for docking fingerprint computation.

TODO(bramsundar): Most of the files below cannot be meaningfully tested
without a body of PDBs that exhibit the various amino-acids in questions.
Build up such a body of PDBs which we can use in the test work.

The code below is a heavily modified version of from Jacob Durrant's
NNScore 2.0.1. The following notice is copied from the original NNScore
file:
# NNScore 2.01 is released under the GNU General Public License (see
# http://www.gnu.org/licenses/gpl.html).
# If you have any questions, comments, or suggestions, please don't
# hesitate to contact me, Jacob Durrant, at jdurrant [at] ucsd [dot]
# edu. If you use NNScore 2.01 in your work, please cite [REFERENCE
# HERE].
"""

__author__ = "Bharath Ramsundar"
__license__ = "GNU General Public License"

import math
import textwrap
import numpy as np


class AromaticRing():
  """Holds information about an aromatic ring."""
  def __init__(self, center, indices, plane_coeff, radius):
    """
    Initializes an aromatic.

    Parameters
    ----------
    center: float
      Center of the ring.
    indices: list
      List of the atom indices for ring atoms.
    plane_coeff: list
      A list of elements [a, b, c, d] that define a plane by equation
      a x + b y + c z = d.
    radius: float
      Ring radius from center.
    """
    self.center = center
    self.indices = indices
    # a*x + b*y + c*z = dI think that
    self.plane_coeff = plane_coeff
    self.radius = radius


def average_point(points):
  """Returns the point with averaged coordinates of arguments.

  Parameters
  ----------
  points: list
    List of point objects. 
  Returns
  -------
  pavg: Point object
    Has coordinates the arithmetic average of those of p1 and p2.
  """
  coords = np.array([0, 0, 0])
  for point in points:
    coords[0] += point.x
    coords[1] += point.y
    coords[2] += point.z
  if len(points) > 0:
    return Point(coords=coords/len(points))
  else:
    return Point(coords=coords)


class Point:
  """
  Simple implementation for a point in 3-space.
  """
  x=99999.0
  y=99999.0
  z=99999.0

  def __init__(self, x=None, y=None, z=None, coords=None):
    """
    Inputs can be specified either by explicitly providing x, y, z coords
    or by providing a numpy array of length 3.

    Parameters
    ----------
    x: float
      X-coord.
    y: float
      Y-coord.
    z: float
      Z-coord.
    coords: np.ndarray
      Should be of length 3 in format np.array([x, y, z])
    Raises
    ------
    ValueError: If no arguments are provided.
    """
    if x and y and z:
      self.x, self.y, self.z = x, y, z 
    elif coords is not None:  # Implicit eval doesn't work on numpy arrays.
      self.x, self.y, self.z = coords[0], coords[1], coords[2]
    else:
      raise ValueError("Must specify coordinates for Point!")

  # TODO(bramsundar): Should this be __copy__?
  def copy_of(self):
    return Point(coords=np.array([self.x, self.y, self.z]))

  def dist_to(self, apoint):
    return (math.sqrt(math.pow(self.x - apoint.x,2)
                    + math.pow(self.y - apoint.y,2)
                    + math.pow(self.z - apoint.z,2)))

  def magnitude(self):
    return self.dist_to(Point(coords=np.array([0, 0, 0])))

class Atom:
  """
  Implements a container class for atoms. This class contains useful
  annotations about the atom.
  """

  def __init__ (self, atomname="", residue="",
                coordinates=Point(99999, 99999, 99999), element="",
                PDBIndex="", line="", atomtype="",
                indices_of_atoms_connecting=None, charge=0, resid=0,
                chain="", structure="", comment=""):
    """
    Initializes an atom.

    Assumes that atom is loaded from a PDB file.

    Parameters
    ----------
    atomname: string
      Name of atom. Note that atomname is not the same as residue since
      atomnames often have extra annotations (e.g., CG, NZ, etc).
      TODO(bramsundar): Find a reference for all possible atomnames.
    residue: string:
      Name of protein residue this atom belongs to.
    element: string
      Name of atom's element.
    coordinate: point
      A point object (x, y, z are in Angstroms).
    PDBIndex: string
      Index of the atom in source PDB file.
    line: string
      The line in the PDB file which specifies this atom.
    atomtype: string
      Type of atom (TODO(bramsundar): How is this different from atomname)
    IndicesOfAtomConnecting: list
      The indices (in a PDB object) of all atoms bonded to this one.
    charge: float
      Associated electrostatic charge.
    resid: int
      The residue number in the receptor (listing the protein as a chain from
      N-Terminus to C-Terminus). Assumes this is a protein atom.
    chain: string
      Chain identifier for molecule. See PDB spec.
    structure: string
      One of ALPHA, BETA, or OTHER for the type of protein secondary
      structure this atom resides in (assuming this is a receptor atom).
    comment: string
      Either LIGAND or RECEPTOR depending on whether this is a ligand or
      receptor atom.
    """
    self.atomname = atomname
    self.residue = residue
    self.coordinates = coordinates
    self.element = element
    self.PDBIndex = PDBIndex
    self.line = line
    self.atomtype = atomtype
    if indices_of_atoms_connecting is not None:
      self.indices_of_atoms_connecting = indices_of_atoms_connecting
    else:
      self.indices_of_atoms_connecting = []
    self.charge = charge
    self.resid = resid
    self.chain = chain
    self.structure = structure
    self.comment = comment

  def copy_of(self):
    theatom = Atom()
    theatom.atomname = self.atomname
    theatom.residue = self.residue
    theatom.coordinates = self.coordinates.copy_of()
    theatom.element = self.element
    theatom.PDBIndex = self.PDBIndex
    theatom.line= self.line
    theatom.atomtype= self.atomtype
    theatom.indices_of_atoms_connecting = self.indices_of_atoms_connecting[:]
    theatom.charge = self.charge
    theatom.resid = self.resid
    theatom.chain = self.chain
    theatom.structure = self.structure
    theatom.comment = self.comment

    return theatom

  def create_PDB_line(self, index):
    """
    Generates appropriate ATOM line for pdb file.

    Parameters
    ----------
    index: int
      Index in associated PDB file.
    """
    output = "ATOM "
    output = output + str(index).rjust(6) + self.atomname.rjust(5) + self.residue.rjust(4)
    output = output + ("%.3f" % self.coordinates.x).rjust(18)
    output = output + ("%.3f" % self.coordinates.y).rjust(8)
    output = output + ("%.3f" % self.coordinates.z).rjust(8)
    output = output + self.element.rjust(24)
    return output

  def number_of_neighbors(self):
    """Reports number of neighboring atoms."""
    return len(self.indices_of_atoms_connecting)

  def add_neighbor_atom_indices(self, indices):
    """
    Adds atoms with provided PDB indices as neighbors.

    Parameters
    ----------
    index: list
      List of indices of neighbors in PDB object.
    """
    for index in indices:
      if not (index in self.indices_of_atoms_connecting):
        self.indices_of_atoms_connecting.append(index)

  def side_chain_or_backbone(self): # only really applies to proteins, assuming standard atom names
    if (self.atomname.strip() == "CA" or self.atomname.strip() == "C"
      or self.atomname.strip() == "O" or self.atomname.strip() == "N"):
      return "BACKBONE"
    else:
      return "SIDECHAIN"

  def read_atom_PDB_line(self, Line):
    """
    TODO(rbharath): This method probably belongs in the PDB class, and not
    in the Atom class.

    Reads an ATOM or HETATM line from PDB and instantiates fields.

    Atoms in PDBs are represented by ATOM or HETATM statements. ATOM and
    HETATM statements follow the following record format:

    (see ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf)

    COLUMNS   DATA TYPE       FIELD             DEFINITION
    -------------------------------------------------------------------------------------
    1 - 6     Record name     "ATOM "/"HETATM"
    7 - 11    Integer         serial            Atom serial number.
    13 - 16   Atom            name              Atom name.
    17        Character       altLoc            Alternate location indicator.
    18 - 20   Residue name    resName           Residue name.
    22        Character       chainID           Chain identifier.
    23 - 26   Integer         resSeq            Residue sequence number.
    27        AChar           iCode             Code for insertion of residues.
    31 - 38   Real(8.3)       x                 Orthogonal coordinates for X in Angstroms.
    39 - 46   Real(8.3)       y                 Orthogonal coordinates for Y in Angstroms.
    47 - 54   Real(8.3)       z                 Orthogonal coordinates for Z in Angstroms.
    55 - 60   Real(6.2)       occupancy         Occupancy.
    61 - 66   Real(6.2)       tempFactor        Temperature factor.
    77 - 78   LString(2)      element           Element symbol, right-justified.
    79 - 80   LString(2)      charge            Charge on the atom.
    """
    self.line = Line
    self.atomname = Line[11:16].strip()

    if len(self.atomname)==1:
      self.atomname = self.atomname + "  "
    elif len(self.atomname)==2:
      self.atomname = self.atomname + " "
    elif len(self.atomname)==3:
      # This line is necessary for babel to work, though many PDBs in
      # the PDB would have this line commented out
      self.atomname = self.atomname + " "

    self.coordinates = Point(coords=np.array([float(Line[30:38]),
        float(Line[38:46]), float(Line[46:54])]))

    # now atom type (for pdbqt)
    self.atomtype = Line[77:79].strip().upper()

    if Line[69:76].strip() != "":
      self.charge = float(Line[69:76])
    else:
      self.charge = 0.0

    if self.element == "": # try to guess at element from name
      two_letters = self.atomname[0:2].strip().upper()
      valid_two_letters = ["BR", "CL", "BI", "AS", "AG", "LI",
          "HG", "MG", "MN", "RH", "ZN", "FE"]
      if two_letters in valid_two_letters:
        self.element = two_letters
      else: #So, just assume it's the first letter.
        # Any number needs to be removed from the element name
        self.element = self.atomname
        self.element = self.element.replace('0','')
        self.element = self.element.replace('1','')
        self.element = self.element.replace('2','')
        self.element = self.element.replace('3','')
        self.element = self.element.replace('4','')
        self.element = self.element.replace('5','')
        self.element = self.element.replace('6','')
        self.element = self.element.replace('7','')
        self.element = self.element.replace('8','')
        self.element = self.element.replace('9','')
        self.element = self.element.replace('@','')

        self.element = self.element[0:1].strip().upper()

    self.PDBIndex = Line[6:12].strip()
    self.residue = Line[16:20]
    # this only uses the rightmost three characters, essentially
    # removing unique rotamer identification
    self.residue = " " + self.residue[-3:]

    if Line[23:26].strip() != "": self.resid = int(Line[23:26])
    else: self.resid = 1

    self.chain = Line[21:22]
    if self.residue.strip() == "": self.residue = " MOL"

class Charged():
  """
  A class that represeents a charged atom.
  """
  def __init__(self, coordinates, indices, positive):
    """
    Parameters
    ----------
    coordinates: point
      Coordinates of atom.
    indices: list
      Contains boolean true or false entries for self and neighbors to
      specify if positive or negative charge
    positive: bool
      Whether this atom is positive or negative.
    """
    self.coordinates = coordinates
    self.indices = indices
    self.positive = positive

class PDB:
  """
  PDB file handler class.

  Provides functionality for loading PDB files. Performs a number of
  clean-up and annotation steps (filling in missing bonds, identifying
  aromatic rings, charged groups, and protein secondary structure
  assignation).
  """

  def __init__(self):
    self.all_atoms={}
    self.non_protein_atoms = {}
    self.max_x = -9999.99
    self.min_x = 9999.99
    self.max_y = -9999.99
    self.min_y = 9999.99
    self.max_z = -9999.99
    self.min_z = 9999.99
    self.rotatable_bonds_count = 0
    self.functions = MathFunctions()
    self.protein_resnames = ["ALA", "ARG", "ASN", "ASP", "ASH", "ASX",
      "CYS", "CYM", "CYX", "GLN", "GLU", "GLH", "GLX", "GLY", "HIS",
      "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE",
      "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    self.aromatic_rings = []
    self.charges = [] # a list of objects of type charge (defined below)

  def load_PDB_from_file(self, filename, line_header="", impute_bonds=False):
    """
    Given a PDB file as a list of lines, loads fields into self.

    Parameters
    ----------
    lines: list
      List of PDB file lines.
    line_header: string
      The string header for PDB lines.
    impute_bonds: bool
      Attempt to guess missing bonds between ligand atoms.
    """
    # Reset internal state
    self.__init__()
    # Now load the file into a list
    with open(filename,"r") as f:
      lines = f.readlines()
    self.load_atoms_from_PDB_list(lines)
    self.load_bonds_from_PDB_list(lines)
    self.check_protein_format()
    if impute_bonds:
      self.create_non_protein_atom_bonds_by_distance()
    self.assign_non_protein_aromatic_rings()
    self.assign_protein_aromatic_rings()
    self.assign_non_protein_charges()
    self.assign_protein_charges()


  def load_atoms_from_PDB_list(self, lines):
    """
    Loads atoms from a list of PDB file lines.

    Parameters
    ----------
    lines: list
      List of strings, one per line in original PDB file.
    """
    autoindex = 1
    # going to keep track of atomname_resid_chain pairs, to make sure
    # redundants aren't loaded.  This basically gets rid of rotomers,
    # I think.
    atom_already_loaded = []

    for line in lines:
      if "between atoms" in line and " A " in line:
        self.rotatable_bonds_count = self.rotatable_bonds_count + 1

      if len(line) >= 7:
        # Load atom data (coordinates, etc.)
        if line[0:4]=="ATOM" or line[0:6]=="HETATM":
          cur_atom = Atom()
          cur_atom.read_atom_PDB_line(line)

          # this string unique identifies each atom
          key = (cur_atom.atomname.strip() + "_" +
            str(cur_atom.resid) + "_" + cur_atom.residue.strip() +
            "_" + cur_atom.chain.strip())
          # so this is a receptor atom that has already been loaded once
          if (key in atom_already_loaded
            and cur_atom.residue.strip() in self.protein_resnames):
            print (line_header
                + "WARNING: Duplicate receptor atom detected: \""
                + cur_atom.line.strip() + "\". Not loading this duplicate.")

          # so either the atom hasn't been loaded, or else it's a non-receptor
          # atom so note that non-receptor atoms can have redundant names, but
          # receptor atoms cannot.  This is because protein residues often
          # contain rotamers
          if (not key in atom_already_loaded
              or not cur_atom.residue.strip() in self.protein_resnames):
            # so each atom can only be loaded once. No rotamers.
            atom_already_loaded.append(key)
            # So you're actually reindexing everything here.
            self.all_atoms[autoindex] = cur_atom
            if (not cur_atom.residue[-3:] in self.protein_resnames):
              self.non_protein_atoms[autoindex] = cur_atom

            autoindex = autoindex + 1

  def load_bonds_from_PDB_list(self, lines):
    """
    Loads bonds from list of PDB file lines.

    Bonds in PDBs are represented by CONECT statements. These lines follow
    the following record format:

    (see ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf)

    Columns    DataType    Definition
    ---------------------------------
    1  - 6      String          -
    7  - 11     Int         Atom index.
    12 - 16     Int         Index of bonded atom.
    17 - 21     Int         Index of bonded atom.
    22 - 26     Int         Index of bonded atom.
    27 - 31     Int         Index of bonded atom.

    If more than 4 bonded atoms are present, then a second CONECT record
    must be specified.

    Parameters
    ----------
    lines: List
      List of strings, one per line in original PDB file.
    Raises
    ------
    ValueError: On improperly formatted input.
    """
    for line in lines:
      if "CONECT" in line:
        if len(line) < 31:
          raise ValueError("Bad PDB! "
              "Improperly formatted CONECT line (too short)")
        atom_index = int(line[6:11].strip())
        if atom_index not in self.all_atoms:
          raise ValueError("Bad PDB! "
              "Improper CONECT line: (atom index not loaded)")
        bonded_atoms = []
        ranges = [(11,16), (16,21), (21,26), (26,31)]
        for (lower, upper) in ranges:
          # Check that the range is nonempty.
          if line[lower:upper].strip():
            index = int(line[lower:upper])
            if index not in self.all_atoms:
              raise ValueError("Bad PDB! "
                  "Improper CONECT line: (bonded atom not loaded)")
            bonded_atoms.append(index)
        atom = self.all_atoms[atom_index]
        atom.add_neighbor_atom_indices(bonded_atoms)

  def save_PDB(self, filename):
    """
    Writes a PDB file version of self to filename.

    Parameters
    ----------
    filename: string
      path to desired PDB file output.
    """
    f = open(filename, 'w')
    towrite = self.save_PDB_string()
    # just so no PDB is empty, VMD will load them all
    if towrite.strip() == "":
      towrite = "ATOM      1  X   XXX             0.000   0.000   0.000                       X"
    f.write(towrite)
    f.close()

  def save_PDB_string(self):
    """
    Generates a PDB string version of self. Used by SavePDB.
    """
    ToOutput = ""
    # write coordinates
    for atomindex in self.all_atoms:
      ToOutput = ToOutput + self.all_atoms[atomindex].create_PDB_line(atomindex) + "\n"
    return ToOutput

  def add_new_atom(self, atom):
    """
    Adds an extra atom to this PDB.

    Parameters
    ----------
    atom: object of atom class
      Will be added to self.
    """
    # first get available index
    t = len(self.all_atoms.keys()) + 1

    # now add atom
    self.all_atoms[t] = atom

  def add_new_atoms(self, atoms):
    """
    Convenience function to add many atoms.

    Parameters
    ----------
    atoms: list
      Entries in atoms should be objects of type atom.
    """
    for atom_obj in atoms:
      self.add_new_atom(atom_obj)

  def add_new_non_protein_atom(self, atom):
    """
    Adds an extra non-protein atom to this PDB.

    Parameters
    ----------
    atom: object of atom class
      Will be added to self.
    """
    # first get available index
    t = len(self.all_atoms.keys()) + 1
    # now add atom
    self.all_atoms[t] = atom
    # Add to non-protein list
    self.non_protein_atoms[t] = atom


  def connected_atoms_of_given_element(self, index, con_element):
    """
    Returns indices of all neighbors of atom at index of given elt.

    Parameters
    ----------
    index: integer
      Index of base atom.
    con_element: string
      Name of desired element.
    """
    atom = self.all_atoms[index]
    connected_atoms = []
    for con_index in atom.indices_of_atoms_connecting:
      con_atom = self.all_atoms[con_index]
      if con_atom.element == con_element:
        connected_atoms.append(con_index)
    return connected_atoms

  def connected_heavy_atoms(self, index):
    """
    Returns indices of all connected heavy atoms.

    Parameters
    ----------
    index: integer
      Index of base atom.
    """
    atom = self.all_atoms[index]
    connected_atoms = []
    for con_index in atom.indices_of_atoms_connecting:
      con_atom = self.all_atoms[con_index]
      if con_atom.element != "H":
        connected_atoms.append(con_index)
    return connected_atoms

  def check_protein_format(self):
    """Check that loaded protein structure is self-consistent.

    Helper function called when loading PDB from file.
    """
    curr_res = ""
    first = True
    residue = []

    for atom_index in self.all_atoms:
      atom = self.all_atoms[atom_index]

      key = atom.residue + "_" + str(atom.resid) + "_" + atom.chain

      if first == True:
        curr_res = key
        first = False

      if key != curr_res:

        self.check_protein_format_process_residue(residue, last_key)

        residue = []
        curr_res = key

      residue.append(atom.atomname.strip())
      last_key = key

    self.check_protein_format_process_residue(residue, last_key)

  def print_warning(self, atom, residue, need):
    """
    Prints warning if residue has improper structure.

    Parameters
    ----------
    atom: string
      Name of affected atom.
    residue: string
      Name of affected residue.
    need: string
      Description of need for this atom in residue.
    """
    text = ('WARNING: There is no atom named "%s"' % atom
        + 'in the protein residue ' + last_key + '.'
        + ' Please use standard naming conventions for all'
        + ' protein residues. This atom is needed to determine'
        + ' %s. If this residue is far from the' % need
        + ' active site, this warning may not affect the NNScore.')
    lines = textwrap.wrap(text, 80)
    for line in lines:
      print line
    print

  def check_protein_format_process_residue(self, residue, last_key):
    """
    Check that specified residue in PDB is formatted correctly.

    Parameters
    ----------
    residue: list
      List of atom names in residue.
    last_key: string
      Should be in format RESNAME_RESNUMBER_CHAIN
    """

    resname, resid, chain = last_key.strip().split("_")
    real_resname = resname[-3:]

    if real_resname in self.protein_resnames: # so it's a protein residue

      if not "N" in residue:
        self.print_warning("N", last_key, "secondary structure")
      if not "C" in residue:
        self.print_warning("C", last_key, "secondary structure")
      if not "CA" in residue:
        self.print_warning("CA", last_key, "secondary structure")

      if real_resname == "GLU" or real_resname == "GLH" or real_resname == "GLX":
        if not "OE1" in residue:
          self.print_warning("OE1", last_key, "salt-bridge interactions")
        if not "OE2" in residue:
          self.print_warning("OE2", last_key, "salt-bridge interactions")

      if real_resname == "ASP" or real_resname == "ASH" or real_resname == "ASX":
        if not "OD1" in residue:
          self.print_warning("OD1", last_key, "salt-bridge interactions")
        if not "OD2" in residue:
          self.print_warning("OD2", last_key, "salt-bridge interactions")

      if real_resname == "LYS" or real_resname == "LYN":
        if not "NZ" in residue:
          self.print_warning("NZ", last_key, "pi-cation and salt-bridge interactions")

      if real_resname == "ARG":
        if not "NH1" in residue:
          self.print_warning("NH1", last_key, "pi-cation and salt-bridge interactions")
        if not "NH2" in residue:
          self.print_warning("NH2", last_key, "pi-cation and salt-bridge interactions")

      if real_resname == "HIS" or real_resname == "HID" or real_resname == "HIE" or real_resname == "HIP":
        if not "NE2" in residue:
          self.print_warning("NE2", last_key, "pi-cation and salt-bridge interactions")
        if not "ND1" in residue:
          self.print_warning("ND1", last_key, "pi-cation and salt-bridge interactions")

      if real_resname == "PHE":
        if not "CG" in residue:
          self.print_warning("CG", last_key, "pi-pi and pi-cation interactions")
        if not "CD1" in residue:
          self.print_warning("CD1", last_key, "pi-pi and pi-cation interactions")
        if not "CD2" in residue:
          self.print_warning("CD2", last_key, "pi-pi and pi-cation interactions")
        if not "CE1" in residue:
          self.print_warning("CE1", last_key, "pi-pi and pi-cation interactions")
        if not "CE2" in residue:
          self.print_warning("CE2", last_key, "pi-pi and pi-cation interactions")
        if not "CZ" in residue:
          self.print_warning("CZ", last_key, "pi-pi and pi-cation interactions")

      if real_resname == "TYR":
        if not "CG" in residue:
          self.print_warning("CG", last_key, "pi-pi and pi-cation interactions")
        if not "CD1" in residue:
          self.print_warning("CD1", last_key, "pi-pi and pi-cation interactions")
        if not "CD2" in residue:
          self.print_warning("CD2", last_key, "pi-pi and pi-cation interactions")
        if not "CE1" in residue:
          self.print_warning("CE1", last_key, "pi-pi and pi-cation interactions")
        if not "CE2" in residue:
          self.print_warning("CE2", last_key, "pi-pi and pi-cation interactions")
        if not "CZ" in residue:
          self.print_warning("CZ", last_key, "pi-pi and pi-cation interactions")

      if real_resname == "TRP":
        if not "CG" in residue:
          self.print_warning("CG", last_key, "pi-pi and pi-cation interactions")
        if not "CD1" in residue:
          self.print_warning("CD1", last_key, "pi-pi and pi-cation interactions")
        if not "CD2" in residue:
          self.print_warning("CD2", last_key, "pi-pi and pi-cation interactions")
        if not "NE1" in residue:
          self.print_warning("NE1", last_key, "pi-pi and pi-cation interactions")
        if not "CE2" in residue:
          self.print_warning("CE2", last_key, "pi-pi and pi-cation interactions")
        if not "CE3" in residue:
          self.print_warning("CE3", last_key, "pi-pi and pi-cation interactions")
        if not "CZ2" in residue:
          self.print_warning("CZ2", last_key, "pi-pi and pi-cation interactions")
        if not "CZ3" in residue:
          self.print_warning("CZ3", last_key, "pi-pi and pi-cation interactions")
        if not "CH2" in residue:
          self.print_warning("CH2", last_key, "pi-pi and pi-cation interactions")

      if (real_resname == "HIS" or real_resname == "HID" or
        real_resname == "HIE" or real_resname == "HIP"):
        if not "CG" in residue:
          self.print_warning("CG", last_key, "pi-pi and pi-cation interactions")
        if not "ND1" in residue:
          self.print_warning("ND1", last_key, "pi-pi and pi-cation interactions")
        if not "CD2" in residue:
          self.print_warning("CD2", last_key, "pi-pi and pi-cation interactions")
        if not "CE1" in residue:
          self.print_warning("CE2", last_key, "pi-pi and pi-cation interactions")
        if not "NE2" in residue:
          self.print_warning("NE2", last_key, "pi-pi and pi-cation interactions")


  # Functions to determine the bond connectivity based on distance
  # ==============================================================

  def create_non_protein_atom_bonds_by_distance(self):
    """
    Creates bonds between non-protein atoms close to each other in PDB.

    This function is only approximate! Avoid using if possible.
    """
    for AtomIndex1 in self.non_protein_atoms:
      atom1 = self.non_protein_atoms[AtomIndex1]
      if not atom1.residue[-3:] in self.protein_resnames: # so it's not a protein residue
        for AtomIndex2 in self.non_protein_atoms:
          if AtomIndex1 != AtomIndex2:
            atom2 = self.non_protein_atoms[AtomIndex2]
            if not atom2.residue[-3:] in self.protein_resnames: # so it's not a protein residue
              dist = self.functions.distance(atom1.coordinates, atom2.coordinates)
              if (dist < self.bond_length(atom1.element, atom2.element) * 1.2):
                print "Adding bond between %d and %d" % (AtomIndex1, AtomIndex2)
                print "Distance between %d and %d is %f" % (AtomIndex1, AtomIndex2, dist)
                atom1.add_neighbor_atom_indices([AtomIndex2])
                atom2.add_neighbor_atom_indices([AtomIndex1])

  def bond_length(self, element1, element2):
    """
    Returns approximate bond-length between atoms of element1 and element2.

    Bond lengths taken from Handbook of Chemistry and Physics.  The
    information provided there was very specific, so representative
    examples were used to specify the bond lengths.  Sitautions could
    arise where these lengths would be incorrect, probably slight errors
    (<0.06) in the hundreds.

    Parameters
    ----------
    element1: string:
      Name of first element.
    element2: string
      Name of second element.
    """
    # All distances are in Angstroms. Duplicate pairs not specified. For
    # example, to find distance ("H", "C"), the lookup key is ("C", "H")
    distances = {
      ("C", "C"): 1.53,
      ("N", "N"): 1.425,
      ("O", "O"): 1.469,
      ("S", "S"): 2.048,
      ("SI", "SI"): 2.359,

      ("C", "H"): 1.059,
      ("C", "N"): 1.469,
      ("C", "O"): 1.413,
      ("C", "S"): 1.819,
      ("C", "F"): 1.399,
      ("C", "CL"): 1.790,
      ("C", "BR"): 1.910,
      ("C", "I"): 2.162,

      ("N", "H"): 1.009,
      ("N", "O"): 1.463,
      ("N", "BR"): 1.843,
      ("N", "CL"): 1.743,
      ("N", "F"): 1.406,
      ("N", "I"): 2.2,

      ("O", "S"): 1.577,
      ("O", "H"): 0.967,

      # This one not from source sited above. Not sure where it's from, but
      # it wouldn't ever be used in the current context ("AutoGrow")
      ("S", "H"): 2.025/1.5,
      ("S", "N"): 1.633,
      ("S", "BR"): 2.321,
      ("S", "CL"): 2.283,
      ("S", "F"): 1.640,
      ("S", "I"): 2.687,

      ("P", "BR"): 2.366,
      ("P", "CL"): 2.008,
      ("P", "F"): 1.495,
      ("P", "I"): 2.490,
      # estimate based on eye balling Handbook of Chemistry and Physics
      ("P", "O"): 1.6,


      ("SI", "BR"): 2.284,
      ("SI", "CL"): 2.072,
      ("SI", "F"): 1.636,
      ("SI", "P"): 2.264,
      ("SI", "S"): 2.145,
      ("SI", "C"): 1.888,
      ("SI", "N"): 1.743,
      ("SI", "O"): 1.631,

      ("H", "H"): .7414,
    }
    if (element1, element2) in distances:
      return distances[(element1, element2)]
    elif (element2, element1) in distances:
      return distances[(element2, element1)]
    else:
      raise ValueError("Distance between %s and %s is unknown" %
          (element1, element2))

  # Functions to identify positive charges
  # ======================================

  def identify_metallic_charges(self):
    """Assign charges to metallic ions.

    Returns
    -------
    charges: list
      Contains a Charge object for every metallic cation.
    """
    # Metallic atoms are assumed to be cations.
    charges = []
    for atom_index in self.non_protein_atoms:
      atom = self.non_protein_atoms[atom_index]
      if (atom.element == "MG" or atom.element == "MN" or
          atom.element == "RH" or atom.element == "ZN" or
          atom.element == "FE" or atom.element == "BI" or
          atom.element == "AS" or atom.element == "AG"):
        chrg = Charged(atom.coordinates, [atom_index], True)
        charges.append(chrg)
    return charges

  def identify_nitrogen_group_charges(self):
    """Assign charges to nitrogen groups where necessary.

    Returns
    -------
    charges: list
      Contains a Charge object for every charged nitrogen group.
    """
    charges = []
    for atom_index in self.non_protein_atoms:
      atom = self.non_protein_atoms[atom_index]
      # Get all the quartenary amines on non-protein residues (these are the
      # only non-protein groups that will be identified as positively
      # charged). Note that nitrogen has only 5 valence electrons (out of 8
      # for a full shell), so any nitrogen with four bonds must be positively
      # charged (think NH4+).
      if atom.element == "N":
        # a quartenary amine, so it's easy
        if atom.number_of_neighbors() == 4:
          indexes = [atom_index]
          indexes.extend(atom.indices_of_atoms_connecting)
          # so the indices stored is just the index of the nitrogen and any
          # attached atoms
          chrg = Charged(atom.coordinates, indexes, True)
          charges.append(chrg)
        # maybe you only have two hydrogens added, but they're sp3 hybridized.
        # Just count this as a quartenary amine, since I think the positive
        # charge would be stabilized. This situation can arise with
        # lone-pair electron nitrogen compounds like pyrrolidine
        # (http://www.chem.ucla.edu/harding/tutorials/lone_pair.pdf)
        elif atom.number_of_neighbors() == 3:
          nitrogen = atom
          atom1 = self.all_atoms[atom.indices_of_atoms_connecting[0]]
          atom2 = self.all_atoms[atom.indices_of_atoms_connecting[1]]
          atom3 = self.all_atoms[atom.indices_of_atoms_connecting[2]]
          angle1 = (self.functions.angle_between_three_points(atom1.coordinates,
            nitrogen.coordinates, atom2.coordinates) * 180.0 / math.pi)
          angle2 = (self.functions.angle_between_three_points(atom1.coordinates,
            nitrogen.coordinates, atom3.coordinates) * 180.0 / math.pi)
          angle3 = (self.functions.angle_between_three_points(atom2.coordinates,
            nitrogen.coordinates, atom3.coordinates) * 180.0 / math.pi)
          average_angle = (angle1 + angle2 + angle3) / 3
          # Test that the angles approximately match the tetrahedral 109
          # degrees
          if math.fabs(average_angle - 109.0) < 5.0:
            indexes = [atom_index]
            indexes.extend(atom.indices_of_atoms_connecting)
            # so indexes added are the nitrogen and any attached atoms.
            chrg = Charged(nitrogen.coordinates, indexes, True)
            charges.append(chrg)
    return charges

  def identify_phosphorus_group_charges(self):
    """Assign charges to phosphorus groups where necessary.

    Searches for phosphate-like groups and assigns charges.
    
    Returns
    -------
    charges: list
      Contains a Charge object for every charged phosphorus group.
    """
    charges = []
    for atom_index in self.non_protein_atoms:
      atom = self.non_protein_atoms[atom_index]
      # let's check for a phosphate or anything where a phosphorus is bound
      # to two oxygens, where both oxygens are bound to only one heavy atom
      # (the phosphorus). I think this will get several phosphorus
      # substances.
      if atom.element == "P":
        oxygens = self.connected_atoms_of_given_element(atom_index,"O")
        if len(oxygens) >=2: # the phosphorus is bound to at least two oxygens
          # now count the number of oxygens that are only bound to the phosphorus
          count = 0
          for oxygen_index in oxygens:
            if len(self.connected_heavy_atoms(oxygen_index)) == 1: count = count + 1
          if count >=2: # so there are at least two oxygens that are only bound to the phosphorus
            indexes = [atom_index]
            indexes.extend(oxygens)
            chrg = Charged(atom.coordinates, indexes, False)
            charges.append(chrg)
    return charges

  def identify_carbon_group_charges(self):
    """Assign charges to carbon groups where necessary.

    Checks for guanidino-like groups and carboxylates.

    TODO(rbharath): This function is monolithic and very special-purpose.
    Can some more general design be created here?

    Returns
    -------
    charges: list
      Contains a Charge object for every charged carbon group.
    """
    charges = []
    for atom_index in self.non_protein_atoms:
      atom = self.non_protein_atoms[atom_index]
      # let's check for guanidino-like groups (actually H2N-C-NH2,
      # where not CN3.)
      if atom.element == "C":
        # if the carbon has only three atoms connected to it
        if atom.number_of_neighbors() == 3:
          nitrogens = self.connected_atoms_of_given_element(atom_index, "N")
          # if true, carbon is connected to at least two nitrogens now,
          # so we need to count the number of nitrogens that are only
          # connected to one heavy atom (the carbon)
          if len(nitrogens) >= 2:
            print "Found carbon with two nitrogens attached!"

            nitrogens_to_use = []
            all_connected = atom.indices_of_atoms_connecting[:]
            # Index of atom that connects this charged group to
            # the rest of the molecule, ultimately to make sure
            # it's sp3 hybridized. Remains -1 if no such atom exists.
            connector_ind = -1

            for atmindex in nitrogens:
              if len(self.connected_heavy_atoms(atmindex)) == 1:
                nitrogens_to_use.append(atmindex)
                all_connected.remove(atmindex)

            # TODO(rbharath): Is picking the first non-nitrogen atom
            # correct here?
            if len(all_connected) > 0:
              connector_ind = all_connected[0]

            # Handle case of guanidinium cation
            if len(nitrogens_to_use) == 3 and connector_ind == -1:
              pt = atom.coordinates.copy_of()
              charges.append(Charged(pt, [atom_index], True))
            elif len(nitrogens_to_use) == 2 and connector_ind != -1:
              # so there are at two nitrogens that are only
              # connected to the carbon (and probably some
              # hydrogens)

              # now you need to make sure connector_ind atom is sp3 hybridized
              connector_atom = self.all_atoms[connector_ind]
              if ((connector_atom.element == "C" and
                  connector_atom.number_of_neighbors() == 4)
                or (connector_atom.element == "O"
                  and connector_atom.number_of_neighbors() == 2)
                or connector_atom.element == "N"
                or connector_atom.element == "S"
                or connector_atom.element == "P"):

                # There are only two "guanidino" nitrogens. Assume the
                # negative charge is spread equally between the two.
                avg_pt = average_point(
                    [self.all_atoms[nitrogen].coordinates for nitrogen in
                     nitrogens_to_use])

                indexes = [atom_index]
                indexes.extend(nitrogens_to_use)
                indexes.extend(self.connected_atoms_of_given_element(nitrogens_to_use[0],"H"))
                indexes.extend(self.connected_atoms_of_given_element(nitrogens_to_use[1],"H"))

                charges.append(Charged(avg_pt, indexes, True)) # True because it's positive

      if atom.element == "C": # let's check for a carboxylate
          # a carboxylate carbon will have three items connected to it.
          if atom.number_of_neighbors() == 3:
            oxygens = self.connected_atoms_of_given_element(atom_index, "O")
            # a carboxylate will have two oxygens connected to
            # it. Now, each of the oxygens should be connected
            # to only one heavy atom (so if it's connected to a
            # hydrogen, that's okay)
            if len(oxygens) == 2:
              if (len(self.connected_heavy_atoms(oxygens[0])) == 1
                and len(self.connected_heavy_atoms(oxygens[1])) == 1):
                # so it's a carboxylate! Add a negative charge.

                # Assume negative charge is centered between the two
                # oxygens.
                avg_pt = average_point(
                    [self.all_atoms[oxygen].coordinates for oxygen in
                    oxygens])
                chrg = Charged(avg_pt,
                    [oxygens[0], atom_index, oxygens[1]], False)
                charges.append(chrg)
    return charges

  def identify_sulfur_group_charges(self):
    """Assigns charges to sulfur groups.

    Searches for Sulfonates.

    Returns
    -------
    charges: list
      Contains a Charge object for every charged sulfur group.
    """
    charges = []
    for atom_index in self.non_protein_atoms:
      atom = self.non_protein_atoms[atom_index]
      # let's check for a sulfonate or anything where a sulfur is
      # bound to at least three oxygens and at least three are
      # bound to only the sulfur (or the sulfur and a hydrogen).
      if atom.element == "S":
        oxygens = self.connected_atoms_of_given_element(atom_index,"O")
        # the sulfur is bound to at least three oxygens now
        # count the number of oxygens that are only bound to the
        # sulfur
        if len(oxygens) >=3:
          count = 0
          for oxygen_index in oxygens:
            if len(self.connected_heavy_atoms(oxygen_index)) == 1: count = count + 1
          # so there are at least three oxygens that are only
          # bound to the sulfur
          if count >=3:
            indexes = [atom_index]
            indexes.extend(oxygens)
            chrg = Charged(atom.coordinates, indexes, False)
            charges.append(chrg)
    return charges


  def assign_non_protein_charges(self):
    """
    Assign positive and negative charges to non-protein atoms.

    This function handles the following cases:

      1) Metallic ions (assumed to be cations)
      2) Quartenary amines (such as NH4+)
      2) sp3 hybridized nitrogen (such as pyrrolidine)
      3) Carboxylates (RCOO-)
      4) Guanidino Groups (NHC(=NH)NH2)
      5) Phosphates (PO4(3-))
      6) Sulfonate (RSO2O-)
    """
    self.charges += self.identify_metallic_charges()
    self.charges += self.identify_nitrogen_group_charges()
    self.charges += self.identify_carbon_group_charges()
    self.charges += self.identify_phosphorus_group_charges()
    self.charges += self.identify_sulfur_group_charges()

  def get_residues(self):
    """Returns a list of all residues in this protein.

    This function uses keys of the following type to uniquely identify
    protein residues: RESNAME_RESNUMBER_CHAIN.

    Returns
    -------
    residues: list
      Each list element is a tuple whose first element is a key (of type
      defined above) and whose second element is a list of the atom-indices
      that make up this residue.
    """
    keys = []
    cur_key = None
    cur_res = []
    residues = []
    for atom_index in self.all_atoms:
      atom = self.all_atoms[atom_index]
      # Assign each atom a residue key.
      key = atom.residue + "_" + str(atom.resid) + "_" + atom.chain
      if not cur_key:
        cur_key = key

      if key != cur_key:
        keys.append(cur_key)
        residues.append(cur_res)
        cur_key = key
        cur_res = []

      cur_res.append(atom_index)
    # Handle edge case of last residue.
    keys.append(cur_key)
    residues.append(cur_res)
    return zip(keys, residues)
      

  def assign_protein_charges(self):
    """Assigns charges to atoms in charged residues.

    This function uses keys of the following type to uniquely identify
    protein residues: RESNAME_RESNUMBER_CHAIN.
    """
    res_list = self.get_residues()
    self.charges += self.get_lysine_charges(res_list)
    self.charges += self.get_arginine_charges(res_list)
    self.charges += self.get_histidine_charges(res_list)
    self.charges += self.get_glutamic_acid_charges(res_list)
    self.charges += self.get_aspartic_acid_charges(res_list)

  def get_residue_charges(self, res_list, resnames, atomnames,
      charged_atomnames, positive=True):
    """Assign charges to specified residue.

    Regardless of protonation state, we assume below that residues are
    charged, since evidence in the literature ("The Cation Pi Interaction,"
    TODO(rbharath): Verify citation) suggests that charges will be
    stabilized.

    Parameters
    ---------
    res_list: list
      List of tuples output by get_residue_list
    resnames: list
      List of acceptable names for residue (e.g. [PHE], [HIS, HIP, HIE,
      HID])
    atomnames: list
      List of names of atoms in charged group.
    charged_atomnames: list
      List of atoms which will be averaged to yield charge location.
    positive: bool
      Whether charge is positive or not.
    Returns
    -------
    aromatics: list
      List of Aromatic objects.
    """
    charges = []
    num_matches = 0
    for key, res in res_list:
      resname, resid, chain = key.strip().split("_")
      real_resname = resname[-3:]
      if real_resname in resnames:
        num_matches += 1
        indices = []
        charged_atoms = [] # The terminal nitrogen holds charge.
        for index in res:
          atom = self.all_atoms[index]
          atomname = atom.atomname.strip()
          if atomname in atomnames:
            indices.append(index)
          if atomname in charged_atomnames:
            charged_atoms.append(atom)
        if len(charged_atoms) == len(charged_atomnames):
          avg_pt = average_point([n.coordinates for n in
              charged_atoms])
          if avg_pt.magnitude() != 0:
            charges.append(Charged(avg_pt, indices, positive))
    print "NUM_MATCHES = %d" % num_matches
    return charges

  def get_lysine_charges(self, res_list):
    """Assign charges to lysine residues.
    
    Regardless of protonation state, assume that lysine is charged.
    Recall that LYS is positive charged lysine and LYN is neutral. See
    http://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/addh/addh.html

    TODO(rbharath): The get_*_charge functions are all highly similar.
    Refactor out the common logic and make each of these functions an
    invocation of that code.

    Parameters
    ----------
    res_list: list
      List of tuples output by get_residue_list
    """
    return self.get_residue_charges(res_list, ["LYS", "LYN"],
        ["NZ", "HZ1", "HNZ1", "HZ2", "HNZ2", "HZ3", "HNZ3"],
        ["NZ"])
    #charges = []
    #for key, res in res_list:
    #  resname, resid, chain = key.strip().split("_")
    #  real_resname = resname[-3:]
    #  if real_resname == "LYS" or real_resname == "LYN":
    #    indices = []
    #    nitrogen_atom = None  # The terminal nitrogen holds charge.
    #    for index in res:
    #      atom = self.all_atoms[index]
    #      atomname = atom.atomname.strip()
    #      if atomname == "NZ":
    #        indices.append(index)
    #        nitrogen_atom = atom
    #      if atomname == "HZ1" or atomname == "HNZ1":
    #        indices.append(index)
    #      if atomname == "HZ2" or atomname == "HNZ2":
    #        indices.append(index)
    #      if atomname == "HZ3" or atomname == "HNZ3":
    #        indices.append(index)
    #    charges.append(Charged(nitrogen_atom.coordinates, indices, True))
    #return charges

  def get_arginine_charges(self, res_list):
    """Assign charges to arginine residues.

    Parameters
    ----------
    res_list: list
      List of tuples output by get_residue_list
    """
    return self.get_residue_charges(res_list, ["ARG"],
        ["NH1", "NH2", "2HH2", "HN22", "1HH2", "HN12", "2HH1", "HN21",
        "1HH1", "HN11"], ["NH1", "NH2"])
    #charges = []
    #for key, res in res_list:
    #  resname, _, _ = key.strip().split("_")
    #  real_resname = resname[-3:]
    #  if real_resname == "ARG":
    #    indices = []
    #    charged_nitrogens = []
    #    for index in res:
    #      atom = self.all_atoms[index]
    #      atomname = atom.atomname.strip()
    #      if atomname == "NH1":
    #        charged_nitrogens.append(atom)
    #        indices.append(index)
    #      if atomname == "NH2":
    #        charged_nitrogens.append(atom)
    #        indices.append(index)
    #      if atomname == "2HH2" or atomname == "HN22":
    #        indices.append(index)
    #      if atomname == "1HH2" or atomname == "HN12":
    #        indices.append(index)
    #      if atomname == "CZ":
    #        indices.append(index)
    #      if atomname == "2HH1" or atomname == "HN21":
    #        indices.append(index)
    #      if atomname == "1HH1" or atomname == "HN11":
    #        indices.append(index)

    #    if len(charged_nitrogens) == 2:
    #      avg_pt = average_point([n.coordinates for n in
    #          charged_nitrogens])
    #      if avg_pt.magnitude() != 0:
    #        charges.append(Charged(avg_pt, indices, True))
    #return charges

  def get_histidine_charges(self, res_list):
    """Assign charges to histidine residues.

    The specific histidine name determines the protonation state:

    * HID: Protonate delta-Nitrogen.
    * HIE: Protonate epsilon-Nitrogen.
    * HIP: Protonate both nitrogens.
    * HIS: Protonation unspecified.

    Regardless of protonation state, assume it's charged. This based on
    "The Cation-Pi Interaction," which suggests protonated state would
    be stabilized. But let's not consider HIS when doing salt bridges.
    Parameters
    ----------
    res_list: list
      List of tuples output by get_residue_list
    """
    return self.get_residue_charges(res_list, ["HIS", "HID", "HIE", "HIP"],
        ["NE2", "ND1", "HE2", "HD1", "CE1", "CD2", "CG"],
        ["NE2", "ND1"])
    #charges = []
    #for key, res in res_list:
    #  resname, _, _ = key.strip().split("_")
    #  real_resname = resname[-3:]
    #  if (real_resname == "HIS" or real_resname == "HID" or
    #    real_resname == "HIE" or real_resname == "HIP"):
    #    indices = []
    #    charged_nitrogens = []
    #    for index in res:
    #      atom = self.all_atoms[index]
    #      atomname = atom.atomname.strip()
    #      if atomname == "NE2":
    #        charged_nitrogens.append(atom)
    #        indices.append(index)
    #      if atomname == "ND1":
    #        charged_nitrogens.append(atom)
    #        indices.append(index)
    #      if (atomname == "HE2" or atomname == "HD1"
    #       or atomname == "CE1" or atomname == "CD2"
    #       or atomname == "CG"):
    #        indices.append(index)

    #    if len(charged_nitrogens) == 2:
    #      avg_pt = average_point([n.coordinates for n in
    #          charged_nitrogens])
    #      if avg_pt.magnitude() != 0:
    #        charges.append(Charged(avg_pt, indices, True))
    #return charges

  def get_glutamic_acid_charges(self, res_list):
    """Assign charges to histidine residues.

    The specific glutamic acid name determines the protonation state:

    * GLU: Negatively charged (deprotonated).
    * GLH: Neutral charge (protonated).
    * GLX: Protonation unspecified.

    See
    http://aria.pasteur.fr/documentation/use-aria/version-2.2/non-standard-atom-or-residue-definitions
    or 
    http://proteopedia.org/wiki/index.php/Standard_Residues

    Regardless of protonation state, assume it's charged. This based on
    "The Cation-Pi Interaction," which suggests protonated state would
    be stabilized..

    Parameters
    ----------
    res_list: list
      List of tuples output by get_residue_list
    """
    return self.get_residue_charges(res_list, ["GLU", "GLH", "GLX"],
        ["OE1", "OE2", "CD"], ["OE1", "OE2"], positive=False)
    #charges = []
    #for key, res in res_list:
    #  resname, _, _ = key.strip().split("_")
    #  real_resname = resname[-3:]
    #  if real_resname == "GLU" or real_resname == "GLH" or real_resname == "GLX":
    #    # regardless of protonation state, assume it's charged. This based on
    #    # "The Cation-Pi Interaction," which suggests protonated state would
    #    # be stabilized.
    #    indices = []
    #    charged_oxygens = []
    #    for index in res:
    #      atom = self.all_atoms[index]
    #      atomname = atom.atomname.strip()
    #      if atomname == "OE1":
    #        charged_oxygens.append(atom)
    #        indices.append(index)
    #      if atomname == "OE2":
    #        charged_oxygens.append(atom)
    #        indices.append(index)
    #      if atomname == "CD":
    #        indices.append(index)

    #    if len(charged_oxygens) == 2:
    #      avg_pt = average_point([n.coordinates for n in charged_oxygens])
    #      if avg_pt.magnitude() != 0:
    #        charges.append(Charged(avg_pt, indices, False))
    #return charges


  def get_aspartic_acid_charges(self, res_list):
    """Assign charges to aspartic acid residues.

    The specific aspartic acid name determines the protonation.

    * ASP: Negatively charged (deprotonated).
    * ASH: Neutral charge (protonated).
    * ASX: Protonation unspecified.

    Regardless of protonation state, assume it's charged. This based on
    "The Cation-Pi Interaction," which suggests protonated state would
    be stabilized.
    Parameters
    ----------
    res_list: list
      List of tuples output by get_residue_list
    """
    return self.get_residue_charges(res_list, ["ASP", "ASH", "ASX"],
        ["OD1", "OD2", "CG"], ["OD1", "OD2"], positive=False)
    #charges = []
    #for key, res in res_list:
    #  resname, resid, chain = key.strip().split("_")
    #  real_resname = resname[-3:]

    #  # TODO(bramsundar): This comment about Cation-Pi interactions
    #  # is repeated in multiple places. Look into this interaction
    #  # and verify that it holds true for the residues in question.
    #  if (real_resname == "ASP" or real_resname == "ASH" or
    #    real_resname == "ASX"):
    #    charge_pt = Point(coords=np.array([0.0,0.0,0.0]))
    #    count = 0.0
    #    indices = []
    #    charged_oxygens = []
    #    for index in res:
    #      atom = self.all_atoms[index]
    #      if atom.atomname.strip() == "OD1":
    #        charged_oxygens.append(atom)
    #        indices.append(index)
    #      if atom.atomname.strip() == "OD2":
    #        charged_oxygens.append(atom)
    #        indices.append(index)
    #      if atom.atomname.strip() == "CG": indices.append(index)

    #    if len(charged_oxygens) == 2:
    #      avg_pt = average_point([n.coordinates for n in charged_oxygens])
    #      if avg_pt.magnitude() != 0:
    #        charges.append(Charged(avg_pt, indices, False))
    #return charges


  # Functions to identify aromatic rings
  # ====================================

  def get_aromatic_marker(self, indices_of_ring):
    """Identify aromatic markers.

    Parameters
    ----------
    indices_of_ring: list
      Contains atom indices for all atoms in the ring.
    """
    # first identify the center point
    points_list = []
    total = len(indices_of_ring)
    x_sum = 0.0
    y_sum = 0.0
    z_sum = 0.0

    for index in indices_of_ring:
      atom = self.all_atoms[index]
      points_list.append(atom.coordinates)
      x_sum = x_sum + atom.coordinates.x
      y_sum = y_sum + atom.coordinates.y
      z_sum = z_sum + atom.coordinates.z

    if total == 0:
      return # to prevent errors in some cases

    center = Point(coords=np.array([x_sum / total, y_sum / total, z_sum /
        total]))

    # now get the radius of the aromatic ring
    radius = 0.0
    for index in indices_of_ring:
      atom = self.all_atoms[index]
      dist = center.dist_to(atom.coordinates)
      if dist > radius:
        radius = dist

    # now get the plane that defines this ring
    if len(indices_of_ring) < 3:
      # to prevent an error in some cases. If there aren't three point, you can't define a plane
      return
    elif len(indices_of_ring) == 3:
      A = self.all_atoms[indices_of_ring[0]].coordinates
      B = self.all_atoms[indices_of_ring[1]].coordinates
      C = self.all_atoms[indices_of_ring[2]].coordinates
    elif len(indices_of_ring) == 4:
      A = self.all_atoms[indices_of_ring[0]].coordinates
      B = self.all_atoms[indices_of_ring[1]].coordinates
      C = self.all_atoms[indices_of_ring[3]].coordinates
    else: # best, for 5 and 6 member rings
      A = self.all_atoms[indices_of_ring[0]].coordinates
      B = self.all_atoms[indices_of_ring[2]].coordinates
      C = self.all_atoms[indices_of_ring[4]].coordinates

    AB = self.functions.vector_subtraction(B,A)
    AC = self.functions.vector_subtraction(C,A)
    ABXAC = self.functions.CrossProduct(AB,AC)

    # formula for plane will be ax + by + cz = d
    x1 = self.all_atoms[indices_of_ring[0]].coordinates.x
    y1 = self.all_atoms[indices_of_ring[0]].coordinates.y
    z1 = self.all_atoms[indices_of_ring[0]].coordinates.z

    a = ABXAC.x
    b = ABXAC.y
    c = ABXAC.z
    d = a*x1 + b*y1 + c*z1

    ar_ring = AromaticRing(center, indices_of_ring, [a,b,c,d], radius)
    #self.aromatic_rings.append(ar_ring)
    return ar_ring

  def assign_non_protein_aromatic_rings(self):
    """Identifies aromatic rings in ligand atoms.

    TODO(rbharath): This function is still monolithic. Better refactoring?
    """
    # Get all the rings containing each of the atoms in the ligand
    all_rings = []
    for atom_index in self.non_protein_atoms:
      all_rings.extend(self.all_rings_containing_atom(atom_index))

    # Ensure that no ring is a subset of another.
    # TODO(rbharath): When is this ever the case?
    for ring_index_1 in range(len(all_rings)):
      ring1 = all_rings[ring_index_1]
      if len(ring1) != 0:
        for ring_index_2 in range(len(all_rings)):
          if ring_index_1 != ring_index_2:
            ring2 = all_rings[ring_index_2]
            if len(ring2) != 0:
              if self.set1_is_subset_of_set2(ring1, ring2) == True:
                all_rings[ring_index_2] = []

    while [] in all_rings:
      all_rings.remove([])

    # Now we need to figure out which of these ligands are aromatic
    # (planar)
    for ring_index in range(len(all_rings)):
      ring = all_rings[ring_index]
      is_flat = True
      for t in range(-3, len(ring)-3):
        pt1 = self.non_protein_atoms[ring[t]].coordinates
        pt2 = self.non_protein_atoms[ring[t+1]].coordinates
        pt3 = self.non_protein_atoms[ring[t+2]].coordinates
        pt4 = self.non_protein_atoms[ring[t+3]].coordinates

        # first, let's see if the last atom in this ring is a carbon
        # connected to four atoms. That would be a quick way of
        # telling this is not an aromatic ring
        cur_atom = self.non_protein_atoms[ring[t+3]]
        if cur_atom.element == "C" and cur_atom.number_of_neighbors() == 4:
          is_flat = False
          break

        # now check the dihedral between the ring atoms to see if
        # it's flat
        angle = self.functions.dihedral(pt1, pt2, pt3, pt4) * 180 / math.pi
        # 15 degrees is the cutoff, ring[t], ring[t+1], ring[t+2],
        # ring[t+3] range of this function is -pi to pi
        if (angle > -165 and angle < -15) or (angle > 15 and angle < 165):
          is_flat = False
          break

        # now check the dihedral between the ring atoms and an atom
        # connected to the current atom to see if that's flat too.
        for substituent_atom_index in cur_atom.indices_of_atoms_connecting:
          pt_sub = self.non_protein_atoms[substituent_atom_index].coordinates
          angle = self.functions.dihedral(pt2, pt3, pt4, pt_sub) * 180 / math.pi
          # 15 degress is the cutoff, ring[t], ring[t+1], ring[t+2],
          # ring[t+3], range of this function is -pi to pi
          if (angle > -165 and angle < -15) or (angle > 15 and angle < 165):
            is_flat = False
            break

      if is_flat == False:
        all_rings[ring_index] = []
      # While I'm at it, three and four member rings are not aromatic
      if len(ring) < 5:
        all_rings[ring_index] = []
      # While I'm at it, if the ring has more than 6, also throw it out. So
      # only 5 and 6 member rings are allowed.
      if len(ring) > 6:
        all_rings[ring_index] = []

    while [] in all_rings:
      all_rings.remove([])

    for ring in all_rings:
      self.aromatic_rings.append(self.get_aromatic_marker(ring))

  def all_rings_containing_atom(self, index):
    """Identify all rings that contain atom at index.

    Parameters
    ----------
    index: int
      Index of provided atom.
    """

    all_rings = []

    atom = self.all_atoms[index]
    for connected_atom in atom.indices_of_atoms_connecting:
      self.ring_recursive(connected_atom, [index], index, all_rings)

    return all_rings

  def ring_recursive(self, index, already_crossed, orig_atom, all_rings):
    """Recursive helper function for ring identification.

    Parameters
    ----------
    index: int
      Index of specified atom.
    already_crossed: list
      TODO(rbharath)
    orig_atom: int
      Index of the original atom in ring.
    all_rings: list
      Used to recursively build up ring structure.
    """

    if len(already_crossed) > 6:
      # since you're only considering aromatic rings containing 5 or 6
      # members anyway, save yourself some time.
      return

    atom = self.all_atoms[index]

    temp = already_crossed[:]
    temp.append(index)

    for connected_atom in atom.indices_of_atoms_connecting:
      if not connected_atom in already_crossed:
        self.ring_recursive(connected_atom, temp, orig_atom, all_rings)
      if connected_atom == orig_atom and orig_atom != already_crossed[-1]:
        all_rings.append(temp)


  def assign_protein_aromatic_rings(self):
    """Identifies aromatic rings in protein residues.
    """
    res_list = self.get_residues()
    self.aromatic_rings += self.get_phenylalanine_aromatics(res_list)
    self.aromatic_rings += self.get_tyrosine_aromatics(res_list)
    self.aromatic_rings += self.get_histidine_aromatics(res_list)
    self.aromatic_rings += self.get_tryptophan_aromatics(res_list)
    #for key, res in res_list:
    #  self.assign_aromatic_rings_from_protein_process_residue(res, key)

  def get_residue_aromatics(self, res_list, resname, ring_atomnames):
    """Helper function that identifies aromatics in given residue.

    Parameters
    ----------
    res_list: list
      List of tuples output by get_residue_list
    resname: list
      List of acceptable names for residue (e.g. [PHE], [HIS, HIP, HIE,
      HID])
    ring_atomnames: list
      List of names of atoms in aromatic ring.
    Returns
    -------
    aromatics: list
      List of Aromatic objects.
    """
    aromatics = []
    num_matches = 0
    print "resname: %s" % resname
    for key, res in res_list:
      real_resname, resid, chain = key.strip().split("_")[-3:]
      indices_of_ring = []
      if real_resname in resname:
        num_matches += 1
        indices_of_ring = []
        for index in res:
          if self.all_atoms[index].atomname.strip() in ring_atomnames:
            indices_of_ring.append(index)
        aromatics.append(self.get_aromatic_marker(indices_of_ring))
    print "num_matches: %d" % num_matches
    return aromatics


  def get_phenylalanine_aromatics(self, res_list):
    """Assign aromatics in phenylalanines.
    
    Parameters
    ----------
    res_list: list
      List of tuples output by get_residue_list
    Returns
    -------
    aromatics: list
      List of Aromatic objects for aromatics in phenylalanines.
    """
    return self.get_residue_aromatics(res_list, "PHE",
        ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"])

  def get_tyrosine_aromatics(self, res_list):
    """Assign aromatics in tyrosines.
    
    Parameters
    ----------
    res_list: list
      List of tuples output by get_residue_list
    Returns
    -------
    aromatics: list
      List of Aromatic objects for aromatics in tyrosines.
    """
    return self.get_residue_aromatics(res_list, "TYR",
        ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"])

  def get_histidine_aromatics(self, res_list):
    """Assign aromatics in histidines.
    
    Parameters
    ----------
    res_list: list
      List of tuples output by get_residue_list
    Returns
    -------
    aromatics: list
      List of Aromatic objects for aromatics in histidines.
    """
    return self.get_residue_aromatics(res_list,
        ["HIS", "HID", "HIE", "HIP"],
        ["CG", "ND1", "CE1", "NE2", "CD2"])

  def get_tryptophan_aromatics(self, res_list):
    """Assign aromatics in tryptophans.
    
    Parameters
    ----------
    res_list: list
      List of tuples output by get_residue_list
    Returns
    -------
    aromatics: list
      List of Aromatic objects for aromatics in tryptophans.
    """
    # Tryptophan has two aromatic rings.
    small_ring = self.get_residue_aromatics(res_list,
        ["TRP"],
        ["CG", "CD1", "NE1", "CE2", "CD2"])
    large_ring = self.get_residue_aromatics(res_list,
        ["TRP"],
        ["CE2", "CD2", "CE3", "CZ3", "CH2", "CZ2"])
    return small_ring + large_ring


  #def assign_aromatic_rings_from_protein_process_residue(self, residue, last_key):
  #  """Find aromatic rings in protein residues.

  #  Parameters
  #  ---------
  #  residue: list
  #    List of atom indices for this residue.
  #  last_key: string
  #    keys have format RESNAME_RESNUMBER_CHAIN
  #  """
  #  resname, resid, chain = last_key.strip().split("_")
  #  real_resname = resname[-3:]

  #  #if real_resname == "PHE":

  #  #  # Note that order is important in the following.
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CG": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CD1": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CE1": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CZ": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CE2": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CD2": indices_of_ring.append(index)

  #  #  self.aromatic_rings.append(self.get_aromatic_marker(indices_of_ring))

  #  #if real_resname == "TYR":
  #  #  indices_of_ring = []

  #  #  # Note that order is important in the following.
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CG": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CD1": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CE1": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CZ": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CE2": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CD2": indices_of_ring.append(index)

  #  #  self.aromatic_rings.append(self.get_aromatic_marker(indices_of_ring))

  #  #if (real_resname == "HIS"
  #  #  or real_resname == "HID"
  #  #  or real_resname == "HIE"
  #  #  or real_resname == "HIP"):
  #  #  indices_of_ring = []

  #  #  # Note that order is important in the following.
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CG": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "ND1": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CE1": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "NE2": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CD2": indices_of_ring.append(index)

  #  #  self.aromatic_rings.append(self.get_aromatic_marker(indices_of_ring))

  #  #if real_resname == "TRP":
  #  #  indices_of_ring = []

  #  #  # Note that order is important in the following.
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CG": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CD1": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "NE1": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CE2": indices_of_ring.append(index)
  #  #  for index in residue:
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CD2": indices_of_ring.append(index)

  #  #  self.aromatic_rings.append(self.get_aromatic_marker(indices_of_ring))

  #  #  indices_of_ring = []

  #  #  for index in residue: # written this way because order is important
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CE2": indices_of_ring.append(index)
  #  #  for index in residue: # written this way because order is important
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CD2": indices_of_ring.append(index)
  #  #  for index in residue: # written this way because order is important
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CE3": indices_of_ring.append(index)
  #  #  for index in residue: # written this way because order is important
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CZ3": indices_of_ring.append(index)
  #  #  for index in residue: # written this way because order is important
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CH2": indices_of_ring.append(index)
  #  #  for index in residue: # written this way because order is important
  #  #    atom = self.all_atoms[index]
  #  #    if atom.atomname.strip() == "CZ2": indices_of_ring.append(index)

  #  #  self.aromatic_rings.append(self.get_aromatic_marker(indices_of_ring))

  # TODO(bramsundar): This looks like it should be a standard
  # python function.
  def set1_is_subset_of_set2(self, set1, set2):
    is_subset = True
    for item in set1:
      if not item in set2:
        is_subset = False
        break
    return is_subset

  # Functions to assign secondary structure to protein residues
  # ===========================================================

  def assign_secondary_structure(self):
    """Assign secondary structure labels (assuming self is a protein).

    keys in this function have form RESNUMBER_CHAIN where CHAIN is
    the chain identifier for this molecule.

    TODO(bramsundar): This function is monolithic. Can we break it down?
    """
    # first, we need to know what resid's are available
    resids = []
    last_key = "-99999_Z"
    for atom_index in self.all_atoms:
      atom = self.all_atoms[atom_index]
      key = str(atom.resid) + "_" + atom.chain
      if key != last_key:
        last_key = key
        resids.append(last_key)

    structure = {}
    for resid in resids:
      structure[resid] = "OTHER"

    atoms = []

    for atom_index in self.all_atoms:
      atom = self.all_atoms[atom_index]
      if atom.side_chain_or_backbone() == "BACKBONE":
        if len(atoms) < 8:
          atoms.append(atom)
        else:
          atoms.pop(0)
          atoms.append(atom)

          # now make sure the first four all have the same resid and
          # the last four all have the same resid
          if (atoms[0].resid == atoms[1].resid
            and atoms[0].resid == atoms[2].resid
            and atoms[0].resid == atoms[3].resid
            and atoms[0] != atoms[4].resid
            and atoms[4].resid == atoms[5].resid
            and atoms[4].resid == atoms[6].resid
            and atoms[4].resid == atoms[7].resid
            and atoms[0].resid + 1 == atoms[7].resid
            and atoms[0].chain == atoms[7].chain):

            resid1 = atoms[0].resid
            resid2 = atoms[7].resid

            # Now give easier to use names to the atoms
            for atom in atoms:
              if atom.resid == resid1 and atom.atomname.strip() == "N":
                first_N = atom
              if atom.resid == resid1 and atom.atomname.strip() == "C":
                first_C = atom
              if atom.resid == resid1 and atom.atomname.strip() == "CA":
                first_CA = atom

              if atom.resid == resid2 and atom.atomname.strip() == "N":
                second_N = atom
              if atom.resid == resid2 and atom.atomname.strip() == "C":
                second_C = atom
              if atom.resid == resid2 and atom.atomname.strip() == "CA":
                second_CA = atom

            # Now compute the phi and psi dihedral angles
            phi = self.functions.dihedral(first_C.coordinates, second_N.coordinates,
                second_CA.coordinates, second_C.coordinates) * 180.0 / math.pi
            psi = self.functions.dihedral(first_N.coordinates, first_CA.coordinates,
                first_C.coordinates, second_N.coordinates) * 180.0 / math.pi

            # Now use those angles to determine if it's alpha or beta
            if phi > -145 and phi < -35 and psi > -70 and psi < 50:
              key1 = str(first_C.resid) + "_" + first_C.chain
              key2 = str(second_C.resid) + "_" + second_C.chain
              structure[key1] = "ALPHA"
              structure[key2] = "ALPHA"
            # beta. This gets some loops (by my eye), but it's the best I could do.
            if ((phi >= -180 and phi < -40 and psi <= 180 and psi > 90)
              or (phi >= -180 and phi < -70 and psi <= -165)):
              key1 = str(first_C.resid) + "_" + first_C.chain
              key2 = str(second_C.resid) + "_" + second_C.chain
              structure[key1] = "BETA"
              structure[key2] = "BETA"

    # Now update each of the atoms with this structural information
    for atom_index in self.all_atoms:
      atom = self.all_atoms[atom_index]
      key = str(atom.resid) + "_" + atom.chain
      atom.structure = structure[key]

    # Some more post processing.
    CA_list = [] # first build a list of the indices of all the alpha carbons
    for atom_index in self.all_atoms:
      atom = self.all_atoms[atom_index]
      if (atom.residue.strip() in self.protein_resnames
        and atom.atomname.strip() == "CA"):
        CA_list.append(atom_index)

    # some more post processing.
    change = True
    while change == True:
      change = False

      # A residue of index i is only going to be in an alpha helix
      # its CA is within 6 A of the CA of the residue i + 3
      for CA_atom_index in CA_list:
        CA_atom = self.all_atoms[CA_atom_index]
        if CA_atom.structure == "ALPHA":
          # so it's in an alpha helix
          another_alpha_is_close = False
          for other_CA_atom_index in CA_list:
            # so now compare that CA to all the other CA's
            other_CA_atom = self.all_atoms[other_CA_atom_index]
            if other_CA_atom.structure == "ALPHA": # so it's also in an alpha helix
              if other_CA_atom.resid - 3 == CA_atom.resid or other_CA_atom.resid + 3 == CA_atom.resid:
                # so this CA atom is one of the ones the first atom
                # might hydrogen bond with
                if other_CA_atom.coordinates.dist_to(CA_atom.coordinates) < 6.0:
                  # so these two CA atoms are close enough together
                  # that their residues are probably hydrogen bonded
                  another_alpha_is_close = True
                  break
          if another_alpha_is_close == False:
            self.set_structure_of_residue(CA_atom.chain, CA_atom.resid, "OTHER")
            change = True

      # Alpha helices are only alpha helices if they span at least 4
      # residues (to wrap around and hydrogen bond). I'm going to
      # require them to span at least 5 residues, based on
      # examination of many structures.
      for index_in_list in range(len(CA_list)-5):

        index_in_pdb1 = CA_list[index_in_list]
        index_in_pdb2 = CA_list[index_in_list+1]
        index_in_pdb3 = CA_list[index_in_list+2]
        index_in_pdb4 = CA_list[index_in_list+3]
        index_in_pdb5 = CA_list[index_in_list+4]
        index_in_pdb6 = CA_list[index_in_list+5]

        atom1 = self.all_atoms[index_in_pdb1]
        atom2 = self.all_atoms[index_in_pdb2]
        atom3 = self.all_atoms[index_in_pdb3]
        atom4 = self.all_atoms[index_in_pdb4]
        atom5 = self.all_atoms[index_in_pdb5]
        atom6 = self.all_atoms[index_in_pdb6]

        if (atom1.resid + 1 == atom2.resid
          and atom2.resid + 1 == atom3.resid
          and atom3.resid + 1 == atom4.resid
          and atom4.resid + 1 == atom5.resid
          and atom5.resid + 1 == atom6.resid): # so they are sequential
          if (atom1.structure != "ALPHA"
            and atom2.structure == "ALPHA"
            and atom3.structure != "ALPHA"):
            self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
            change = True
          if (atom2.structure != "ALPHA"
            and atom3.structure == "ALPHA"
            and atom4.structure != "ALPHA"):
            self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
            change = True
          if (atom3.structure != "ALPHA"
            and atom4.structure == "ALPHA"
            and atom5.structure != "ALPHA"):
            self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
            change = True
          if (atom4.structure != "ALPHA"
            and atom5.structure == "ALPHA"
            and atom6.structure != "ALPHA"):
            self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
            change = True

          if (atom1.structure != "ALPHA"
            and atom2.structure == "ALPHA"
            and atom3.structure == "ALPHA"
            and atom4.structure != "ALPHA"):
            self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
            self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
            change = True
          if (atom2.structure != "ALPHA"
            and atom3.structure == "ALPHA"
            and atom4.structure == "ALPHA"
            and atom5.structure != "ALPHA"):
            self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
            self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
            change = True
          if (atom3.structure != "ALPHA"
            and atom4.structure == "ALPHA"
            and atom5.structure == "ALPHA"
            and atom6.structure != "ALPHA"):
            self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
            self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
            change = True

          if (atom1.structure != "ALPHA"
            and atom2.structure == "ALPHA"
            and atom3.structure == "ALPHA"
            and atom4.structure == "ALPHA"
            and atom5.structure != "ALPHA"):
            self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
            self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
            self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
            change = True
          if (atom2.structure != "ALPHA"
            and atom3.structure == "ALPHA"
            and atom4.structure == "ALPHA"
            and atom5.structure == "ALPHA"
            and atom6.structure != "ALPHA"):
            self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
            self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
            self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
            change = True

          if (atom1.structure != "ALPHA"
            and atom2.structure == "ALPHA"
            and atom3.structure == "ALPHA"
            and atom4.structure == "ALPHA"
            and atom5.structure == "ALPHA"
            and atom6.structure != "ALPHA"):
            self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
            self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
            self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
            self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
            change = True

      # now go through each of the BETA CA atoms. A residue is only
      # going to be called a beta sheet if CA atom is within 6.0 A
      # of another CA beta, same chain, but index difference > 2.
      for CA_atom_index in CA_list:
          CA_atom = self.all_atoms[CA_atom_index]
          if CA_atom.structure == "BETA":
            # so it's in a beta sheet
            another_beta_is_close = False
            for other_CA_atom_index in CA_list:
              if other_CA_atom_index != CA_atom_index:
                # so not comparing an atom to itself
                other_CA_atom = self.all_atoms[other_CA_atom_index]
                if other_CA_atom.structure == "BETA":
                  # so you're comparing it only to other BETA-sheet atoms
                  if other_CA_atom.chain == CA_atom.chain:
                    # so require them to be on the same chain. needed to indecies can be fairly compared
                    if math.fabs(other_CA_atom.resid - CA_atom.resid) > 2:
                      # so the two residues are not simply adjacent to each other on the chain
                      if CA_atom.coordinates.dist_to(other_CA_atom.coordinates) < 6.0:
                        # so these to atoms are close to each other
                        another_beta_is_close = True
                        break
            if another_beta_is_close == False:
              self.set_structure_of_residue(CA_atom.chain, CA_atom.resid, "OTHER")
              change = True

      # Now some more post-processing needs to be done. Do this
      # again to clear up mess that may have just been created
      # (single residue beta strand, for example)
      # Beta sheets are usually at least 3 residues long

      for index_in_list in range(len(CA_list)-3):

        index_in_pdb1 = CA_list[index_in_list]
        index_in_pdb2 = CA_list[index_in_list+1]
        index_in_pdb3 = CA_list[index_in_list+2]
        index_in_pdb4 = CA_list[index_in_list+3]

        atom1 = self.all_atoms[index_in_pdb1]
        atom2 = self.all_atoms[index_in_pdb2]
        atom3 = self.all_atoms[index_in_pdb3]
        atom4 = self.all_atoms[index_in_pdb4]

        if (atom1.resid + 1 == atom2.resid and atom2.resid + 1 ==
          atom3.resid and atom3.resid + 1 == atom4.resid):
          # so they are sequential

          if (atom1.structure != "BETA"
            and atom2.structure == "BETA"
            and atom3.structure != "BETA"):
            self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
            change = True
          if (atom2.structure != "BETA"
            and atom3.structure == "BETA"
            and atom4.structure != "BETA"):
            self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
            change = True
          if (atom1.structure != "BETA"
            and atom2.structure == "BETA"
            and atom3.structure == "BETA"
            and atom4.structure != "BETA"):
            self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
            self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
            change = True

  def set_structure_of_residue(self, chain, resid, structure):
    for atom_index in self.all_atoms:
      atom = self.all_atoms[atom_index]
      if atom.chain == chain and atom.resid == resid:
        atom.structure = structure

# TODO(bramsundar): Rework this to use numpy instead of custom
# implementations.
class MathFunctions:

  def vector_subtraction(self, vector1, vector2): # vector1 - vector2
    return Point(coords=np.array([vector1.x - vector2.x, vector1.y -
        vector2.y, vector1.z - vector2.z]))

  def CrossProduct(self, Pt1, Pt2): # never tested
    Response = Point(coords=np.array([0,0,0]))

    Response.x = Pt1.y * Pt2.z - Pt1.z * Pt2.y
    Response.y = Pt1.z * Pt2.x - Pt1.x * Pt2.z
    Response.z = Pt1.x * Pt2.y - Pt1.y * Pt2.x

    return Response;

  def vector_scalar_multiply(self, vector, scalar):
    return Point(coords=np.array([vector.x * scalar, vector.y * scalar,
    vector.z * scalar]))

  def dot_product(self, point1, point2):
    return point1.x * point2.x + point1.y * point2.y + point1.z * point2.z

  def dihedral(self, point1, point2, point3, point4): # never tested
    """Compute dihedral angle between 4 points."""

    b1 = self.vector_subtraction(point2, point1)
    b2 = self.vector_subtraction(point3, point2)
    b3 = self.vector_subtraction(point4, point3)

    b2Xb3 = self.CrossProduct(b2,b3)
    b1Xb2 = self.CrossProduct(b1,b2)

    b1XMagb2 = self.vector_scalar_multiply(b1,b2.magnitude())
    radians = math.atan2(self.dot_product(b1XMagb2,b2Xb3), self.dot_product(b1Xb2,b2Xb3))
    return radians

  def angle_between_three_points(self, point1, point2, point3): # As in three connected atoms
    vector1 = self.vector_subtraction(point1, point2)
    vector2 = self.vector_subtraction(point3, point2)
    return self.angle_between_points(vector1, vector2)

  def angle_between_points(self, point1, point2):
    new_point1 = self.return_normalized_vector(point1)
    new_point2 = self.return_normalized_vector(point2)
    dot_prod = self.dot_product(new_point1, new_point2)
    if dot_prod > 1.0: dot_prod = 1.0 # to prevent errors that can rarely occur
    if dot_prod < -1.0: dot_prod = -1.0
    return math.acos(dot_prod)

  def return_normalized_vector(self, vector):
    dist = self.distance(Point(coords=np.array([0,0,0])), vector)
    return Point(coords=np.array([vector.x/dist, vector.y/dist,
        vector.z/dist]))

  def distance(self, point1, point2):
    return point1.dist_to(point2)

  def project_point_onto_plane(self, apoint, plane_coefficients):
    # essentially finds the point on the plane that is closest to the
    # specified point the plane_coefficients are [a,b,c,d], where the
    # plane is ax + by + cz = d

    # First, define a plane using cooeficients a, b, c, d such that ax + by + cz = d
    a = plane_coefficients[0]
    b = plane_coefficients[1]
    c = plane_coefficients[2]
    d = plane_coefficients[3]

    # Now, define a point in space (s,u,v)
    s = apoint.x
    u = apoint.y
    v = apoint.z

    # the formula of a line perpendicular to the plan passing through (s,u,v) is:
    #x = s + at
    #y = u + bt
    #z = v + ct

    t = (d - a*s - b*u - c*v) / (a*a + b*b + c*c)

    # here's the point closest on the plane
    x = s + a*t
    y = u + b*t
    z = v + c*t

    return Point(coords=np.array([x,y,z]))
