"""
Helper Classes and Functions for NNScoreFingerprint Computation.

The code below is adopted from Jacob Durrant's NNScore 2.0.1. The
following notive is copied from the original NNScore file:
# NNScore 2.01 is released under the GNU General Public License (see
# http://www.gnu.org/licenses/gpl.html).
# If you have any questions, comments, or suggestions, please don't
# hesitate to contact me, Jacob Durrant, at jdurrant [at] ucsd [dot]
# edu. If you use NNScore 2.01 in your work, please cite [REFERENCE
# HERE].
"""

# TODO(bramsundar): Should I make these classes have capital names to
# maintain consistency with other source files?

__author__ = "Bharath Ramsundar"
__license__ = "GNU General Public License"

import math
import textwrap


class point:
  """
  Simple implementation for a point in 3-space.

  TODO(bramsundar): Rework this implementation to use numpy vectors.
  """
  x=99999.0
  y=99999.0
  z=99999.0

  def __init__ (self, x, y ,z):
    self.x = x
    self.y = y
    self.z = z

  # TODO(bramsundar): Should this be __copy__?
  def copy_of(self):
    return point(self.x, self.y, self.z)

  def dist_to(self, apoint):
    return (math.sqrt(math.pow(self.x - apoint.x,2)
                    + math.pow(self.y - apoint.y,2)
                    + math.pow(self.z - apoint.z,2)))

  def Magnitude(self):
    return self.dist_to(point(0,0,0))

  def CreatePDBLine(self, index):
    """
    Generates appropriate ATOM line for pdb file.

    TODO(bramsundar): How is this different from CreatePDBLine for the
    atom class?

    Parameters
    ----------
    index: int
      TODO(bramsundar): What is purpose of index?
    """

    output = "ATOM "
    output = output + str(index).rjust(6) + "X".rjust(5) + "XXX".rjust(4)
    output = output + ("%.3f" % self.x).rjust(18)
    output = output + ("%.3f" % self.y).rjust(8)
    output = output + ("%.3f" % self.z).rjust(8)
    output = output + "X".rjust(24)
    return output

class atom:
  """
  Implements a container class for atoms. This class contains useful
  annotations about the atom.
  """

  def __init__ (self, atomname="", residue="",
                coordinates=point(99999, 99999, 99999), element="",
                PDBIndex="", line="", atomtype="",
                IndicesOfAtomsConnecting=None, charge=0, resid=0,
                chain="", structure="", comment=""):
    """
    Initializes an atom.

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
      TODO(bramsundar): What does this do?
    line: string
      The line in the PDB file which specifies this atom.
    atomtype: string
      Type of atom (TODO(bramsundar): How is this different from atomname)
    IndicesOfAtomConnecting: list
      The indices (in a PDB object) of all atoms bonded to this one.
    charge: float
      Associated electrotstatic charge.
    resid: int
      TODO(bramsundar)
    chain: string
      TODO(bramsundar)
    structure: string
      One of ALPHA, BETA, or OTHER for the type of protein secondary
      structure this atom resides in (assuming this is a receptor atom).
    comments: string
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
    if IndicesOfAtomsConnecting is not None:
      self.IndicesOfAtomsConnecting = IndicesOfAtomsConnecting
    else:
      self.IndicesOfAtomsConnecting = []
    self.charge = charge
    self.resid = resid
    self.chain = chain
    self.structure = structure
    self.comment = comment

  def copy_of(self):
    theatom = atom()
    theatom.atomname = self.atomname
    theatom.residue = self.residue
    theatom.coordinates = self.coordinates.copy_of()
    theatom.element = self.element
    theatom.PDBIndex = self.PDBIndex
    theatom.line= self.line
    theatom.atomtype= self.atomtype
    theatom.IndicesOfAtomsConnecting = self.IndicesOfAtomsConnecting[:]
    theatom.charge = self.charge
    theatom.resid = self.resid
    theatom.chain = self.chain
    theatom.structure = self.structure
    theatom.comment = self.comment

    return theatom

  def CreatePDBLine(self, index):
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

  def NumberOfNeighbors(self):
    return len(self.IndicesOfAtomsConnecting)

  def AddNeighborAtomIndex(self, index):
    """
    Adds atom with PDB index provided as neighbor.

    Parameters
    ----------
    index: int
      Index in PDB object.
    """
    if not (index in self.IndicesOfAtomsConnecting):
      self.IndicesOfAtomsConnecting.append(index)

  def SideChainOrBackBone(self): # only really applies to proteins, assuming standard atom names
    if (self.atomname.strip() == "CA" or self.atomname.strip() == "C"
      or self.atomname.strip() == "O" or self.atomname.strip() == "N"):
      return "BACKBONE"
    else:
      return "SIDECHAIN"

  def ReadPDBLine(self, Line):
    """
    Reads an ATOM line from PDB and instantiates fields.
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

    self.coordinates = point(float(Line[30:38]), float(Line[38:46]), float(Line[46:54]))

    # now atom type (for pdbqt)
    self.atomtype = Line[77:79].strip().upper()

    if Line[69:76].strip() != "":
      self.charge = float(Line[69:76])
    else:
      self.charge = 0.0

    if self.element == "": # try to guess at element from name
      two_letters = self.atomname[0:2].strip().upper()
      if two_letters=='BR':
        self.element='BR'
      elif two_letters=='CL':
        self.element='CL'
      elif two_letters=='BI':
        self.element='BI'
      elif two_letters=='AS':
        self.element='AS'
      elif two_letters=='AG':
        self.element='AG'
      elif two_letters=='LI':
        self.element='LI'
      elif two_letters=='HG':
        self.element='HG'
      elif two_letters=='MG':
        self.element='MG'
      elif two_letters=='MN':
        self.element='MN'
      elif two_letters=='RH':
        self.element='RH'
      elif two_letters=='ZN':
        self.element='ZN'
      elif two_letters=='FE':
        self.element='FE'
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

class charged():
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

## TODO(bramsundar): Uncomment this class once tests are written for
## classes above.
# TODO(bramsundar): Other packages (mdtraj in particular) already
# implement PDB handling. Can this class be reconstructed as a thin
# wrapper around MDTraj or RDKit PDB handling?
class PDB:
  """
  PDB file handler class.

  TODO(bramsundar): Add a a short summary of actual PDB handling
  methodology here.
  """

  def __init__(self):
    self.AllAtoms={}
    self.NonProteinAtoms = {}
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
    self.OrigFileName = ""

  def LoadPDB_from_file(self, FileName, line_header=""):
    """
    Sets fields of self by reading PDB file at path FileName.
    """

    self.line_header=line_header

    # Now load the file into a list
    file = open(FileName,"r")
    lines = file.readlines()
    file.close()
    self.LoadPDB_from_list(lines, self.line_header)

  def LoadPDB_from_list(self, lines, line_header=""):
    """
    Given a PDB file as a list of lines, loads fields into self.

    TODO(bramsundar): Should this be a private method?
    """

    self.line_header=line_header
    autoindex = 1
    self.__init__()
    # going to keep track of atomname_resid_chain pairs, to make sure
    # redundants aren't loaded.  This basically gets rid of rotomers,
    # I think.
    atom_already_loaded = []

    for t in range(0, len(lines)):
      line=lines[t]

      if "between atoms" in line and " A " in line:
        self.rotatable_bonds_count = self.rotatable_bonds_count + 1

      if len(line) >= 7:
        # Load atom data (coordinates, etc.)
        if line[0:4]=="ATOM" or line[0:6]=="HETATM":
          TempAtom = atom()
          TempAtom.ReadPDBLine(line)

          # this string unique identifies each atom
          key = (TempAtom.atomname.strip() + "_" +
            str(TempAtom.resid) + "_" + TempAtom.residue.strip() +
            "_" + TempAtom.chain.strip())
          # so this is a receptor atom that has already been loaded once
          if (key in atom_already_loaded
            and TempAtom.residue.strip() in self.protein_resnames):
            print (self.line_header
                + "WARNING: Duplicate receptor atom detected: \""
                + TempAtom.line.strip() + "\". Not loading this duplicate.")

          # so either the atom hasn't been loaded, or else it's a non-receptor
          # atom so note that non-receptor atoms can have redundant names, but
          # receptor atoms cannot.  This is because protein residues often
          # contain rotamers
          if (not key in atom_already_loaded 
              or not TempAtom.residue.strip() in self.protein_resnames):
            # so each atom can only be loaded once. No rotamers.
            atom_already_loaded.append(key)
            # So you're actually reindexing everything here.
            self.AllAtoms[autoindex] = TempAtom
            if (not TempAtom.residue[-3:] in self.protein_resnames):
              self.NonProteinAtoms[autoindex] = TempAtom

            autoindex = autoindex + 1

    self.CheckProteinFormat()
    # only for the ligand, because bonds can be inferred based on
    # atomnames from PDB
    self.CreateNonProteinAtomBondsByDistance()
    self.assign_aromatic_rings()
    self.AssignNonProteinCharges()
    self.AssignProteinCharges()

  def SavePDB(self, filename):
    """
    Writes a PDB file version of self to filename.

    Parameters
    ----------
    filename: string
      path to desired PDB file output.
    """
    f = open(filename, 'w')
    towrite = self.SavePDBString()
    # just so no PDB is empty, VMD will load them all
    if towrite.strip() == "":
      towrite = "ATOM      1  X   XXX             0.000   0.000   0.000                       X"
    f.write(towrite)
    f.close()

  def SavePDBString(self):
    """
    Generates a PDB string version of self. Used by SavePDB.
    """
    ToOutput = ""
    # write coordinates
    for atomindex in self.AllAtoms:
      ToOutput = ToOutput + self.AllAtoms[atomindex].CreatePDBLine(atomindex) + "\n"
    return ToOutput

  def AddNewAtom(self, atom):
    """
    Adds an extra atom to this PDB.

    Parameters
    ----------
    atom: object of atom class
      Will be added to self.
    """
    # first get available index
    t = len(self.AllAtoms.keys()) + 1

    # now add atom
    self.AllAtoms[t] = atom

  def AddNewAtoms(self, atoms):
    """
    Convenience function to add many atoms.

    Parameters
    ----------
    atoms: list
      Entries in atoms should be objects of type atom.
    """
    for atom_obj in atoms:
      self.AddNewAtom(atom_obj)

  def AddNewNonProteinAtom(self, atom):
    """
    Adds an extra non-protein atom to this PDB.

    Parameters
    ----------
    atom: object of atom class
      Will be added to self.
    """
    # first get available index
    t = len(self.AllAtoms.keys()) + 1
    # now add atom
    self.AllAtoms[t] = atom
    # Add to non-protein list
    self.NonProteinAtoms[t] = atom


  def ConnectedAtomsOfGivenElement(self, index, con_element):
    """
    Returns indices of all neighbors of atom at index of given elt.

    Parameters
    ----------
    index: integer
      Index of base atom.
    con_element: string
      Name of desired element.
    """
    atom = self.AllAtoms[index]
    connected_atoms = []
    for con_index in atom.IndicesOfAtomsConnecting:
      con_atom = self.AllAtoms[con_index]
      if con_atom.element == con_element:
        connected_atoms.append(con_index)
    return connected_atoms

  def ConnectedHeavyAtoms(self, index):
    """
    Returns indices of all connected heavy atoms.

    Parameters
    ----------
    index: integer
      Index of base atom.
    """
    atom = self.AllAtoms[index]
    connected_atoms = []
    for con_index in atom.IndicesOfAtomsConnecting:
      con_atom = self.AllAtoms[con_index]
      if con_atom.element != "H":
        connected_atoms.append(con_index)
    return connected_atoms

#  def CheckProteinFormat(self):
#    curr_res = ""
#    first = True
#    residue = []
#
#    for atom_index in self.AllAtoms:
#      atom = self.AllAtoms[atom_index]
#
#      key = atom.residue + "_" + str(atom.resid) + "_" + atom.chain
#
#      if first == True:
#        curr_res = key
#        first = False
#
#      if key != curr_res:
#
#        self.CheckProteinFormatProcessResidue(residue, last_key)
#
#        residue = []
#        curr_res = key
#
#      residue.append(atom.atomname.strip())
#      last_key = key
#
#    self.CheckProteinFormatProcessResidue(residue, last_key)

  def PrintWarning(self, atom, residue, need):
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

  def CheckProteinFormatProcessResidue(self, residue, last_key):
    """
    Check that specified residue in PDB is formatted correctly.

    Parameters
    ----------
    residue: string 
      Three letter name of residue (i.e., PHE, GLU, etc.) 

    last_key: string
      Should be in format RESIDUENAME_RESIDUENUMBER_CHAIN
    """

    resname, resid, chain = last_key.strip().split("_")
    resname = temp[0]
    real_resname = resname[-3:]
    resid = temp[1]
    chain = temp[2]

    if real_resname in self.protein_resnames: # so it's a protein residue

      if not "N" in residue:
        self.PrintWarning("N", last_key, "secondary structure")
      if not "C" in residue:
        self.PrintWarning("C", last_key, "secondary structure")
      if not "CA" in residue:
        self.PrintWarning("CA", last_key, "secondary structure")

      if real_resname == "GLU" or real_resname == "GLH" or real_resname == "GLX":
        if not "OE1" in residue:
          self.PrintWarning("OE1", last_key, "salt-bridge interactions")
        if not "OE2" in residue:
          self.PrintWarning("OE2", last_key, "salt-bridge interactions")

      if real_resname == "ASP" or real_resname == "ASH" or real_resname == "ASX":
        if not "OD1" in residue:
          self.PrintWarning("OD1", last_key, "salt-bridge interactions")
        if not "OD2" in residue:
          self.PrintWarning("OD2", last_key, "salt-bridge interactions")

      if real_resname == "LYS" or real_resname == "LYN":
        if not "NZ" in residue:
          self.PrintWarning("NZ", last_key, "pi-cation and salt-bridge interactions")

      if real_resname == "ARG":
        if not "NH1" in residue:
          self.PrintWarning("NH1", last_key, "pi-cation and salt-bridge interactions")
        if not "NH2" in residue:
          self.PrintWarning("NH2", last_key, "pi-cation and salt-bridge interactions")

      if real_resname == "HIS" or real_resname == "HID" or real_resname == "HIE" or real_resname == "HIP":
        if not "NE2" in residue:
          self.PrintWarning("NE2", last_key, "pi-cation and salt-bridge interactions")
        if not "ND1" in residue:
          self.PrintWarning("ND1", last_key, "pi-cation and salt-bridge interactions")

      if real_resname == "PHE":
        if not "CG" in residue:
          self.PrintWarning("CG", last_key, "pi-pi and pi-cation interactions")
        if not "CD1" in residue:
          self.PrintWarning("CD1", last_key, "pi-pi and pi-cation interactions")
        if not "CD2" in residue:
          self.PrintWarning("CD2", last_key, "pi-pi and pi-cation interactions")
        if not "CE1" in residue:
          self.PrintWarning("CE1", last_key, "pi-pi and pi-cation interactions")
        if not "CE2" in residue:
          self.PrintWarning("CE2", last_key, "pi-pi and pi-cation interactions")
        if not "CZ" in residue:
          self.PrintWarning("CZ", last_key, "pi-pi and pi-cation interactions")

      if real_resname == "TYR":
        if not "CG" in residue:
          self.PrintWarning("CG", last_key, "pi-pi and pi-cation interactions")
        if not "CD1" in residue:
          self.PrintWarning("CD1", last_key, "pi-pi and pi-cation interactions")
        if not "CD2" in residue:
          self.PrintWarning("CD2", last_key, "pi-pi and pi-cation interactions")
        if not "CE1" in residue:
          self.PrintWarning("CE1", last_key, "pi-pi and pi-cation interactions")
        if not "CE2" in residue:
          self.PrintWarning("CE2", last_key, "pi-pi and pi-cation interactions")
        if not "CZ" in residue:
          self.PrintWarning("CZ", last_key, "pi-pi and pi-cation interactions")

      if real_resname == "TRP":
        if not "CG" in residue:
          self.PrintWarning("CG", last_key, "pi-pi and pi-cation interactions")
        if not "CD1" in residue:
          self.PrintWarning("CD1", last_key, "pi-pi and pi-cation interactions")
        if not "CD2" in residue:
          self.PrintWarning("CD2", last_key, "pi-pi and pi-cation interactions")
        if not "NE1" in residue:
          self.PrintWarning("NE1", last_key, "pi-pi and pi-cation interactions")
        if not "CE2" in residue:
          self.PrintWarning("CE2", last_key, "pi-pi and pi-cation interactions")
        if not "CE3" in residue:
          self.PrintWarning("CE3", last_key, "pi-pi and pi-cation interactions")
        if not "CZ2" in residue:
          self.PrintWarning("CZ2", last_key, "pi-pi and pi-cation interactions")
        if not "CZ3" in residue:
          self.PrintWarning("CZ3", last_key, "pi-pi and pi-cation interactions")
        if not "CH2" in residue:
          self.PrintWarning("CH2", last_key, "pi-pi and pi-cation interactions")

      if (real_resname == "HIS" or real_resname == "HID" or
        real_resname == "HIE" or real_resname == "HIP"):
        if not "CG" in residue:
          self.PrintWarning("CG", last_key, "pi-pi and pi-cation interactions")
        if not "ND1" in residue:
          self.PrintWarning("ND1", last_key, "pi-pi and pi-cation interactions")
        if not "CD2" in residue:
          self.PrintWarning("CD2", last_key, "pi-pi and pi-cation interactions")
        if not "CE1" in residue:
          self.PrintWarning("CE2", last_key, "pi-pi and pi-cation interactions")
        if not "NE2" in residue:
          self.PrintWarning("NE2", last_key, "pi-pi and pi-cation interactions")


  # Functions to determine the bond connectivity based on distance
  # ==============================================================

  def CreateNonProteinAtomBondsByDistance(self):
    """
    Creates bonds between non-protein atoms close to each other in PDB.
    """
    for AtomIndex1 in self.NonProteinAtoms:
      atom1 = self.NonProteinAtoms[AtomIndex1]
      if not atom1.residue[-3:] in self.protein_resnames: # so it's not a protein residue
        for AtomIndex2 in self.NonProteinAtoms:
          if AtomIndex1 != AtomIndex2:
            atom2 = self.NonProteinAtoms[AtomIndex2]
            if not atom2.residue[-3:] in self.protein_resnames: # so it's not a protein residue
              dist = self.functions.distance(atom1.coordinates, atom2.coordinates)

              if (dist < self.BondLength(atom1.element, atom2.element) * 1.2):
                atom1.AddNeighborAtomIndex(AtomIndex2)
                atom2.AddNeighborAtomIndex(AtomIndex1)

  def BondLength(self, element1, element2):
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

  def AssignNonProteinCharges(self):
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
      7) Charged Residues (in helper function)
    """
    AllCharged = []
    for atom_index in self.NonProteinAtoms:
      atom = self.NonProteinAtoms[atom_index]
      # Metallic atoms are assumed to be cations.
      if (atom.element == "MG" or atom.element == "MN" or
          atom.element == "RH" or atom.element == "ZN" or
          atom.element == "FE" or atom.element == "BI" or
          atom.element == "AS" or atom.element == "AG"):
        chrg = charged(atom.coordinates, [atom_index], True)
        self.charges.append(chrg)

      # Get all the quartenary amines on non-protein residues (these are the
      # only non-protein groups that will be identified as positively
      # charged). Note that nitrogen has only 5 valence electrons (out of 8
      # for a full shell), so any nitrogen with four bonds must be positively
      # charged (think NH4+).
      if atom.element == "N":
        # a quartenary amine, so it's easy
        if atom.NumberOfNeighbors() == 4:
          indexes = [atom_index]
          indexes.extend(atom.IndicesOfAtomsConnecting)
          # so the indices stored is just the index of the nitrogen and any attached atoms
          chrg = charged(atom.coordinates, indexes, True)
          self.charges.append(chrg)
        # maybe you only have two hydrogens added, but they're sp3 hybridized.
        # Just count this as a quartenary amine, since I think the positive
        # charge would be stabilized. This situation can arise with
        # lone-pair electron nitrogen compounds like pyrrolidine
        # (http://www.chem.ucla.edu/harding/tutorials/lone_pair.pdf)
        elif atom.NumberOfNeighbors() == 3:
          nitrogen = atom
          atom1 = self.AllAtoms[atom.IndicesOfAtomsConnecting[0]]
          atom2 = self.AllAtoms[atom.IndicesOfAtomsConnecting[1]]
          atom3 = self.AllAtoms[atom.IndicesOfAtomsConnecting[2]]
          angle1 = (self.functions.angle_between_three_points(atom1.coordinates,
            nitrogen.coordinates, atom2.coordinates) * 180.0 /
            math.pi)
          angle2 = (self.functions.angle_between_three_points(atom1.coordinates,
            nitrogen.coordinates, atom3.coordinates) * 180.0 /
            math.pi)
          angle3 = (self.functions.angle_between_three_points(atom2.coordinates,
            nitrogen.coordinates, atom3.coordinates) * 180.0 /
            math.pi)
          average_angle = (angle1 + angle2 + angle3) / 3
          # Test that the angles approximately match the tetrahedral 109
          # degrees
          if math.fabs(average_angle - 109.0) < 5.0:
            indexes = [atom_index]
            indexes.extend(atom.IndicesOfAtomsConnecting)
            # so indexes added are the nitrogen and any attached atoms.
            chrg = charged(nitrogen.coordinates, indexes, True)
            self.charges.append(chrg)

      # let's check for guanidino-like groups (actually H2N-C-NH2,
      # where not CN3.)
      if atom.element == "C":
        # if the carbon has only three atoms connected to it
        if atom.NumberOfNeighbors() == 3:
          nitrogens = self.ConnectedAtomsOfGivenElement(atom_index, "N")
          # if true, carbon is connected to at least two nitrogens now,
          # so we need to count the number of nitrogens that are only
          # connected to one heavy atom (the carbon)
          if len(nitrogens) >= 2:

            nitrogens_to_use = []
            all_connected = atom.IndicesOfAtomsConnecting[:]
            not_isolated = -1

            for atmindex in nitrogens:
              if len(self.ConnectedHeavyAtoms(atmindex)) == 1:
                nitrogens_to_use.append(atmindex)
                all_connected.remove(atmindex)

            if len(all_connected) > 0:
              # get the atom that connects this charged group to
              # the rest of the molecule, ultimately to make sure
              # it's sp3 hybridized
              not_isolated = all_connected[0]

            # so there are at two nitrogens that are only
            # connected to the carbon (and probably some
            # hydrogens)
            if len(nitrogens_to_use) == 2 and not_isolated != -1:

              # now you need to make sure not_isolated atom is sp3 hybridized
              not_isolated_atom = self.AllAtoms[not_isolated]
              if ((not_isolated_atom.element == "C" and
                  not_isolated_atom.NumberOfNeighbors() == 4)
                or (not_isolated_atom.element == "O"
                  and not_isolated_atom.NumberOfNeighbors() == 2)
                or not_isolated_atom.element == "N"
                or not_isolated_atom.element == "S"
                or not_isolated_atom.element == "P"):

                pt = self.AllAtoms[nitrogens_to_use[0]].coordinates.copy_of()
                pt.x = pt.x + self.AllAtoms[nitrogens_to_use[1]].coordinates.x
                pt.y = pt.y + self.AllAtoms[nitrogens_to_use[1]].coordinates.y
                pt.z = pt.z + self.AllAtoms[nitrogens_to_use[1]].coordinates.z
                pt.x = pt.x / 2.0
                pt.y = pt.y / 2.0
                pt.z = pt.z / 2.0

                indexes = [atom_index]
                indexes.extend(nitrogens_to_use)
                indexes.extend(self.ConnectedAtomsOfGivenElement(nitrogens_to_use[0],"H"))
                indexes.extend(self.ConnectedAtomsOfGivenElement(nitrogens_to_use[1],"H"))

                chrg = charged(pt, indexes, True) # True because it's positive
                self.charges.append(chrg)

      if atom.element == "C": # let's check for a carboxylate
          # a carboxylate carbon will have three items connected to it.
          if atom.NumberOfNeighbors() == 3:
            oxygens = self.ConnectedAtomsOfGivenElement(atom_index, "O")
            # a carboxylate will have two oxygens connected to
            # it. Now, each of the oxygens should be connected
            # to only one heavy atom (so if it's connected to a
            # hydrogen, that's okay)
            if len(oxygens) == 2:
              if (len(self.ConnectedHeavyAtoms(oxygens[0])) == 1
                and len(self.ConnectedHeavyAtoms(oxygens[1])) == 1):
                # so it's a carboxylate! Add a negative charge.
                pt = self.AllAtoms[oxygens[0]].coordinates.copy_of()
                pt.x = pt.x + self.AllAtoms[oxygens[1]].coordinates.x
                pt.y = pt.y + self.AllAtoms[oxygens[1]].coordinates.y
                pt.z = pt.z + self.AllAtoms[oxygens[1]].coordinates.z
                pt.x = pt.x / 2.0
                pt.y = pt.y / 2.0
                pt.z = pt.z / 2.0
                chrg = charged(pt, [oxygens[0],
                    atom_index, oxygens[1]], False)
                self.charges.append(chrg)

      # let's check for a phosphate or anything where a phosphorus is bound
      # to two oxygens where both oxygens are bound to only one heavy atom
      # (the phosphorus). I think this will get several phosphorus
      # substances.
      if atom.element == "P":
        oxygens = self.ConnectedAtomsOfGivenElement(atom_index,"O")
        if len(oxygens) >=2: # the phosphorus is bound to at least two oxygens
          # now count the number of oxygens that are only bound to the phosphorus
          count = 0
          for oxygen_index in oxygens:
            if len(self.ConnectedHeavyAtoms(oxygen_index)) == 1: count = count + 1
          if count >=2: # so there are at least two oxygens that are only bound to the phosphorus
            indexes = [atom_index]
            indexes.extend(oxygens)
            chrg = charged(atom.coordinates, indexes, False)
            self.charges.append(chrg)

      # let's check for a sulfonate or anything where a sulfur is
      # bound to at least three oxygens and at least three are
      # bound to only the sulfur (or the sulfur and a hydrogen).
      if atom.element == "S":
        oxygens = self.ConnectedAtomsOfGivenElement(atom_index,"O")
        # the sulfur is bound to at least three oxygens now
        # count the number of oxygens that are only bound to the
        # sulfur
        if len(oxygens) >=3:
          count = 0
          for oxygen_index in oxygens:
            if len(self.ConnectedHeavyAtoms(oxygen_index)) == 1: count = count + 1
          # so there are at least three oxygens that are only
          # bound to the sulfur
          if count >=3:
            indexes = [atom_index]
            indexes.extend(oxygens)
            chrg = charged(atom.coordinates, indexes, False)
            self.charges.append(chrg)


#  def AssignProteinCharges(self):
#    """
#    Assigns charges to atoms in charged residues.
#    """
#    # TODO(bramsundar): Is this comment meant to be here?
#    # Now that you've found all the positive charges in non-protein
#    # residues, it's time to look for aromatic rings in protein
#    # residues
#    curr_res = ""
#    first = True
#    residue = []
#
#    for atom_index in self.AllAtoms:
#      atom = self.AllAtoms[atom_index]
#      key = atom.residue + "_" + str(atom.resid) + "_" + atom.chain
#
#      if first == True:
#        curr_res = key
#        first = False
#
#      if key != curr_res:
#
#        self.AssignChargesFromProteinProcessResidue(residue, last_key)
#        residue = []
#        curr_res = key
#
#      residue.append(atom_index)
#      last_key = key
#
#
#    def AssignChargesFromProteinProcessResidue(self, residue, last_key):
#      temp = last_key.strip().split("_")
#      resname = temp[0]
#      real_resname = resname[-3:]
#      resid = temp[1]
#      chain = temp[2]
#
#      # regardless of protonation state, assume it's charged.
#      if real_resname == "LYS" or real_resname == "LYN":
#        for index in residue:
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "NZ":
#
#            # quickly go through the residue and get the hydrogens
#            # attached to this nitrogen to include in the index list
#            indexes = [index]
#            for index2 in residue:
#                atom2 = self.AllAtoms[index2]
#                if atom2.atomname.strip() == "HZ1": indexes.append(index2)
#                if atom2.atomname.strip() == "HZ2": indexes.append(index2)
#                if atom2.atomname.strip() == "HZ3": indexes.append(index2)
#
#            chrg = charged(atom.coordinates, indexes, True)
#            self.charges.append(chrg)
#        break
#
#      if real_resname == "ARG":
#        charge_pt = point(0.0,0.0,0.0)
#        count = 0.0
#        indices = []
#        for index in residue:
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "NH1":
#            chge_pt.x = charge_pt.x + atom.coordinates.x
#            chge_pt.y = charge_pt.y + atom.coordinates.y
#            chge_pt.z = charge_pt.z + atom.coordinates.z
#            indices.append(index)
#            count = count + 1.0
#          if atom.atomname.strip() == "NH2":
#            chge_pt.x = charge_pt.x + atom.coordinates.x
#            chge_pt.y = charge_pt.y + atom.coordinates.y
#            chge_pt.z = charge_pt.z + atom.coordinates.z
#            indices.append(index)
#            count = count + 1.0
#          if atom.atomname.strip() == "2HH2": indices.append(index)
#          if atom.atomname.strip() == "1HH2": indices.append(index)
#          if atom.atomname.strip() == "CZ": indices.append(index)
#          if atom.atomname.strip() == "2HH1": indices.append(index)
#          if atom.atomname.strip() == "1HH1": indices.append(index)
#
#      # TODO(bramsundar): Formatting here was very messed up. Write a
#      # test to ensure that this function behaves as advertised.
#      if count != 0.0:
#
#        charge_pt.x = charge_pt.x / count
#
#        charge_pt.y = charge_pt.y / count
#        charge_pt.z = charge_pt.z / count
#
#        if charge_pt.x != 0.0 or charge_pt.y != 0.0 or charge_pt.z != 0.0:
#            chrg = charged(charge_pt, indices, True)
#            self.charges.append(chrg)
#
#        if (real_resname == "HIS" or real_resname == "HID" or
#          # regardless of protonation state, assume it's charged. This based on
#          # "The Cation-Pi Interaction," which suggests protonated state would
#          # be stabilized. But let's not consider HIS when doing salt bridges.
#          real_resname == "HIE" or real_resname == "HIP"):
#          charge_pt = point(0.0,0.0,0.0)
#          count = 0.0
#          indices = []
#          for index in residue:
#              atom = self.AllAtoms[index]
#              if atom.atomname.strip() == "NE2":
#                charge_pt.x = charge_pt.x + atom.coordinates.x
#                charge_pt.y = charge_pt.y + atom.coordinates.y
#                charge_pt.z = charge_pt.z + atom.coordinates.z
#                indices.append(index)
#                count = count + 1.0
#      if atom.atomname.strip() == "ND1":
#          charge_pt.x = charge_pt.x + atom.coordinates.x
#          charge_pt.y = charge_pt.y + atom.coordinates.y
#          charge_pt.z = charge_pt.z + atom.coordinates.z
#          indices.append(index)
#          count = count + 1.0
#          if atom.atomname.strip() == "HE2": indices.append(index)
#          if atom.atomname.strip() == "HD1": indices.append(index)
#          if atom.atomname.strip() == "CE1": indices.append(index)
#          if atom.atomname.strip() == "CD2": indices.append(index)
#          if atom.atomname.strip() == "CG": indices.append(index)
#
#          if count != 0.0:
#            charge_pt.x = charge_pt.x / count
#            charge_pt.y = charge_pt.y / count
#            charge_pt.z = charge_pt.z / count
#            if charge_pt.x != 0.0 or charge_pt.y != 0.0 or charge_pt.z != 0.0:
#              chrg = charged(charge_pt, indices, True)
#              self.charges.append(chrg)
#
#        if real_resname == "GLU" or real_resname == "GLH" or real_resname == "GLX":
#          # regardless of protonation state, assume it's charged. This based on
#          # "The Cation-Pi Interaction," which suggests protonated state would
#          # be stabilized.
#          charge_pt = point(0.0,0.0,0.0)
#          count = 0.0
#          indices = []
#          for index in residue:
#              atom = self.AllAtoms[index]
#              if atom.atomname.strip() == "OE1":
#          charge_pt.x = charge_pt.x + atom.coordinates.x
#          charge_pt.y = charge_pt.y + atom.coordinates.y
#          charge_pt.z = charge_pt.z + atom.coordinates.z
#          indices.append(index)
#          count = count + 1.0
#      if atom.atomname.strip() == "OE2":
#        charge_pt.x = charge_pt.x + atom.coordinates.x
#        charge_pt.y = charge_pt.y + atom.coordinates.y
#        charge_pt.z = charge_pt.z + atom.coordinates.z
#        indices.append(index)
#        count = count + 1.0
#        if atom.atomname.strip() == "CD": indices.append(index)
#
#        if count != 0.0:
#          charge_pt.x = charge_pt.x / count
#          charge_pt.y = charge_pt.y / count
#          charge_pt.z = charge_pt.z / count
#          if charge_pt.x != 0.0 or charge_pt.y != 0.0 or charge_pt.z != 0.0:
#            chrg = charged(charge_pt, indices, False) # False because it's a negative charge
#            self.charges.append(chrg)
#
#        # TODO(bramsundar): This comment about Cation-Pi interactions
#        # is repeated in multiple places. Look into this interaction
#        # and verify that it holds true for the residues in question.
#        if (real_resname == "ASP" or real_resname == "ASH" or
#          real_resname == "ASX"):
#          # regardless of protonation state, assume it's charged. This based on
#          # "The Cation-Pi Interaction," which suggests protonated state would
#          # be stabilized.
#          charge_pt = point(0.0,0.0,0.0)
#           count = 0.0
#           indices = []
#           for index in residue:
#             atom = self.AllAtoms[index]
#             if atom.atomname.strip() == "OD1":
#              charge_pt.x = charge_pt.x + atom.coordinates.x
#              charge_pt.y = charge_pt.y + atom.coordinates.y
#              charge_pt.z = charge_pt.z + atom.coordinates.z
#              indices.append(index)
#              count = count + 1.0
#      if atom.atomname.strip() == "OD2":
#        charge_pt.x = charge_pt.x + atom.coordinates.x
#        charge_pt.y = charge_pt.y + atom.coordinates.y
#        charge_pt.z = charge_pt.z + atom.coordinates.z
#        indices.append(index)
#        count = count + 1.0
#        if atom.atomname.strip() == "CG": indices.append(index)
#
#        if count != 0.0:
#            charge_pt.x = charge_pt.x / count
#            charge_pt.y = charge_pt.y / count
#            charge_pt.z = charge_pt.z / count
#            if charge_pt.x != 0.0 or charge_pt.y != 0.0 or charge_pt.z != 0.0:
#              # False because it's a negative charge
#              chrg = charged(charge_pt, indices, False)
#              self.charges.append(chrg)
#
#
#    # Functions to identify aromatic rings
#    # ====================================
#
#    def add_aromatic_marker(self, indices_of_ring):
#        # first identify the center point
#        points_list = []
#        total = len(indices_of_ring)
#        x_sum = 0.0
#        y_sum = 0.0
#        z_sum = 0.0
#
#        for index in indices_of_ring:
#            atom = self.AllAtoms[index]
#            points_list.append(atom.coordinates)
#            x_sum = x_sum + atom.coordinates.x
#            y_sum = y_sum + atom.coordinates.y
#            z_sum = z_sum + atom.coordinates.z
#
#        if total == 0: return # to prevent errors in some cases
#
#        center = point(x_sum / total, y_sum / total, z_sum / total)
#
#        # now get the radius of the aromatic ring
#        radius = 0.0
#        for index in indices_of_ring:
#            atom = self.AllAtoms[index]
#            dist = center.dist_to(atom.coordinates)
#            if dist > radius: radius = dist
#
#        # now get the plane that defines this ring
#        if len(indices_of_ring) < 3:
#          # to prevent an error in some cases. If there aren't three point, you can't define a plane
#          return
#        elif len(indices_of_ring) == 3:
#            A = self.AllAtoms[indices_of_ring[0]].coordinates
#            B = self.AllAtoms[indices_of_ring[1]].coordinates
#            C = self.AllAtoms[indices_of_ring[2]].coordinates
#        elif len(indices_of_ring) == 4:
#            A = self.AllAtoms[indices_of_ring[0]].coordinates
#            B = self.AllAtoms[indices_of_ring[1]].coordinates
#            C = self.AllAtoms[indices_of_ring[3]].coordinates
#        else: # best, for 5 and 6 member rings
#            A = self.AllAtoms[indices_of_ring[0]].coordinates
#            B = self.AllAtoms[indices_of_ring[2]].coordinates
#            C = self.AllAtoms[indices_of_ring[4]].coordinates
#
#        AB = self.functions.vector_subtraction(B,A)
#        AC = self.functions.vector_subtraction(C,A)
#        ABXAC = self.functions.CrossProduct(AB,AC)
#
#        # formula for plane will be ax + by + cz = d
#        x1 = self.AllAtoms[indices_of_ring[0]].coordinates.x
#        y1 = self.AllAtoms[indices_of_ring[0]].coordinates.y
#        z1 = self.AllAtoms[indices_of_ring[0]].coordinates.z
#
#        a = ABXAC.x
#        b = ABXAC.y
#        c = ABXAC.z
#        d = a*x1 + b*y1 + c*z1
#
#        ar_ring = self.aromatic_ring(center, indices_of_ring, [a,b,c,d], radius)
#        self.aromatic_rings.append(ar_ring)
#
#    class aromatic_ring():
#      def __init__(self, center, indices, plane_coeff, radius):
#        self.center = center
#        self.indices = indices
#        # a*x + b*y + c*z = dI think that
#        self.plane_coeff = plane_coeff
#        self.radius = radius
#
#    def assign_aromatic_rings(self):
#      # Get all the rings containing each of the atoms in the ligand
#      AllRings = []
#      for atom_index in self.NonProteinAtoms:
#        AllRings.extend(self.all_rings_containing_atom(atom_index))
#
#      for ring_index_1 in range(len(AllRings)):
#        ring1 = AllRings[ring_index_1]
#        if len(ring1) != 0:
#          for ring_index_2 in range(len(AllRings)):
#            if ring_index_1 != ring_index_2:
#              ring2 = AllRings[ring_index_2]
#              if len(ring2) != 0:
#                if self.set1_is_subset_of_set2(ring1, ring2) == True:
#                  AllRings[ring_index_2] = []
#
#      while [] in AllRings:
#        AllRings.remove([])
#
#      # Now we need to figure out which of these ligands are aromatic
#      # (planar)
#
#      for ring_index in range(len(AllRings)):
#        ring = AllRings[ring_index]
#        is_flat = True
#        for t in range(-3, len(ring)-3):
#          pt1 = self.NonProteinAtoms[ring[t]].coordinates
#          pt2 = self.NonProteinAtoms[ring[t+1]].coordinates
#          pt3 = self.NonProteinAtoms[ring[t+2]].coordinates
#          pt4 = self.NonProteinAtoms[ring[t+3]].coordinates
#
#          # first, let's see if the last atom in this ring is a carbon
#          # connected to four atoms. That would be a quick way of
#          # telling this is not an aromatic ring
#          cur_atom = self.NonProteinAtoms[ring[t+3]]
#          if cur_atom.element == "C" and cur_atom.NumberOfNeighbors() == 4:
#            is_flat = False
#            break
#
#          # now check the dihedral between the ring atoms to see if
#          # it's flat
#          angle = self.functions.dihedral(pt1, pt2, pt3, pt4) * 180 / math.pi
#          # 15 degress is the cutoff, ring[t], ring[t+1], ring[t+2],
#          # ring[t+3] range of this function is -pi to pi
#          if (angle > -165 and angle < -15) or (angle > 15 and angle < 165):
#              is_flat = False
#              break
#
#          # now check the dihedral between the ring atoms and an atom
#          # connected to the current atom to see if that's flat too.
#          for substituent_atom_index in cur_atom.IndicesOfAtomsConnecting:
#              pt_sub = self.NonProteinAtoms[substituent_atom_index].coordinates
#              angle = self.functions.dihedral(pt2, pt3, pt4, pt_sub) * 180 / math.pi
#              # 15 degress is the cutoff, ring[t], ring[t+1], ring[t+2],
#              # ring[t+3], range of this function is -pi to pi
#              if (angle > -165 and angle < -15) or (angle > 15 and angle < 165):
#                  is_flat = False
#                  break
#
#        if is_flat == False: AllRings[ring_index] = []
#        # While I'm at it, three and four member rings are not aromatic
#        if len(ring) < 5: AllRings[ring_index] = []
#        # TODO(bramsundar): What about 7-element rings.
#        # While I'm at it, if the ring has more than 6, also throw it out. So
#        # only 5 and 6 member rings are allowed.
#        if len(ring) > 6: AllRings[ring_index] = []
#
#
#
#      while [] in AllRings: AllRings.remove([])
#
#      for ring in AllRings:
#          self.add_aromatic_marker(ring)
#
#      # Now that you've found all the rings in non-protein residues,
#      # it's time to look for aromatic rings in protein residues
#      curr_res = ""
#      first = True
#      residue = []
#
#      for atom_index in self.AllAtoms:
#        atom = self.AllAtoms[atom_index]
#
#        key = atom.residue + "_" + str(atom.resid) + "_" + atom.chain
#
#        if first == True:
#          curr_res = key
#          first = False
#
#        if key != curr_res:
#
#          self.assign_aromatic_rings_from_protein_process_residue(residue, last_key)
#
#          residue = []
#          curr_res = key
#
#        residue.append(atom_index)
#        last_key = key
#
#      self.assign_aromatic_rings_from_protein_process_residue(residue, last_key)
#
#    def assign_aromatic_rings_from_protein_process_residue(self, residue, last_key):
#      temp = last_key.strip().split("_")
#      resname = temp[0]
#      real_resname = resname[-3:]
#      resid = temp[1]
#      chain = temp[2]
#
#      if real_resname == "PHE":
#        indices_of_ring = []
#
#        # TODO(bramsundar): The comment about ordering is repeated
#        # many times. Factor out.
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CG": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CD1": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CE1": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CZ": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CE2": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CD2": indices_of_ring.append(index)
#
#        self.add_aromatic_marker(indices_of_ring)
#
#      if real_resname == "TYR":
#        indices_of_ring = []
#
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CG": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CD1": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CE1": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CZ": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CE2": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CD2": indices_of_ring.append(index)
#
#        self.add_aromatic_marker(indices_of_ring)
#
#      if real_resname == "HIS" or real_resname == "HID" or real_resname == "HIE" or real_resname == "HIP":
#        indices_of_ring = []
#
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CG": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "ND1": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CE1": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "NE2": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CD2": indices_of_ring.append(index)
#
#        self.add_aromatic_marker(indices_of_ring)
#
#      if real_resname == "TRP":
#        indices_of_ring = []
#
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CG": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CD1": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "NE1": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CE2": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CD2": indices_of_ring.append(index)
#
#        self.add_aromatic_marker(indices_of_ring)
#
#        indices_of_ring = []
#
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CE2": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CD2": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CE3": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CZ3": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CH2": indices_of_ring.append(index)
#        for index in residue: # written this way because order is important
#          atom = self.AllAtoms[index]
#          if atom.atomname.strip() == "CZ2": indices_of_ring.append(index)
#
#        self.add_aromatic_marker(indices_of_ring)
#
#    # TODO(bramsundar): This looks like it should be a standard
#    # python function.
#    def set1_is_subset_of_set2(self, set1, set2):
#      is_subset = True
#      for item in set1:
#        if not item in set2:
#          is_subset = False
#          break
#      return is_subset
#
#    def all_rings_containing_atom(self, index):
#
#      AllRings = []
#
#      atom = self.AllAtoms[index]
#      for conneceted_atom in atom.IndicesOfAtomsConnecting:
#        self.ring_recursive(conneceted_atom, [index], index, AllRings)
#
#      return AllRings
#
#    def ring_recursive(self, index, AlreadyCrossed, orig_atom, AllRings):
#
#      if len(AlreadyCrossed) > 6:
#        # since you're only considering aromatic rings containing 5 or 6
#        # members anyway, save yourself some time.
#        return
#
#      atom = self.AllAtoms[index]
#
#      temp = AlreadyCrossed[:]
#      temp.append(index)
#
#      for conneceted_atom in atom.IndicesOfAtomsConnecting:
#        if not conneceted_atom in AlreadyCrossed:
#          self.ring_recursive(conneceted_atom, temp, orig_atom, AllRings)
#        if conneceted_atom == orig_atom and orig_atom != AlreadyCrossed[-1]:
#          AllRings.append(temp)
#
#    # Functions to assign secondary structure to protein residues
#    # ===========================================================
#
#    def assign_secondary_structure(self):
#      # first, we need to know what resid's are available
#      resids = []
#      last_key = "-99999_Z"
#      for atom_index in self.AllAtoms:
#        atom = self.AllAtoms[atom_index]
#        key = str(atom.resid) + "_" + atom.chain
#        if key != last_key:
#          last_key = key
#          resids.append(last_key)
#
#      structure = {}
#      for resid in resids:
#        structure[resid] = "OTHER"
#
#      atoms = []
#
#      for atom_index in self.AllAtoms:
#        atom = self.AllAtoms[atom_index]
#        if atom.SideChainOrBackBone() == "BACKBONE":
#          if len(atoms) < 8:
#            atoms.append(atom)
#          else:
#            atoms.pop(0)
#            atoms.append(atom)
#
#            # now make sure the first four all have the same resid and
#            # the last four all have the same resid
#            if (atoms[0].resid == atoms[1].resid
#              and atoms[0].resid == atoms[2].resid
#              and atoms[0].resid == atoms[3].resid
#              and atoms[0] != atoms[4].resid
#              and atoms[4].resid == atoms[5].resid
#              and atoms[4].resid == atoms[6].resid
#              and atoms[4].resid == atoms[7].resid
#              and atoms[0].resid + 1 == atoms[7].resid
#              and atoms[0].chain == atoms[7].chain):
#
#              resid1 = atoms[0].resid
#              resid2 = atoms[7].resid
#
#              # Now give easier to use names to the atoms
#              for atom in atoms:
#                if atom.resid == resid1 and atom.atomname.strip() == "N": first_N = atom
#                if atom.resid == resid1 and atom.atomname.strip() == "C": first_C = atom
#                if atom.resid == resid1 and atom.atomname.strip() == "CA": first_CA = atom
#
#                if atom.resid == resid2 and atom.atomname.strip() == "N": second_N = atom
#                if atom.resid == resid2 and atom.atomname.strip() == "C": second_C = atom
#                if atom.resid == resid2 and atom.atomname.strip() == "CA": second_CA = atom
#
#              # Now compute the phi and psi dihedral angles
#              phi = self.functions.dihedral(first_C.coordinates, second_N.coordinates, second_CA.coordinates, second_C.coordinates) * 180.0 / math.pi
#              psi = self.functions.dihedral(first_N.coordinates, first_CA.coordinates, first_C.coordinates, second_N.coordinates) * 180.0 / math.pi
#
#              # Now use those angles to determine if it's alpha or beta
#              if phi > -145 and phi < -35 and psi > -70 and psi < 50:
#                  key1 = str(first_C.resid) + "_" + first_C.chain
#                  key2 = str(second_C.resid) + "_" + second_C.chain
#                  structure[key1] = "ALPHA"
#                  structure[key2] = "ALPHA"
#              # beta. This gets some loops (by my eye), but it's the best I could do.
#              if ((phi >= -180 and phi < -40 and psi <= 180 and psi > 90)
#                or (phi >= -180 and phi < -70 and psi <= -165)):
#                  key1 = str(first_C.resid) + "_" + first_C.chain
#                  key2 = str(second_C.resid) + "_" + second_C.chain
#                  structure[key1] = "BETA"
#                  structure[key2] = "BETA"
#
#      # Now update each of the atoms with this structural information
#      for atom_index in self.AllAtoms:
#        atom = self.AllAtoms[atom_index]
#        key = str(atom.resid) + "_" + atom.chain
#        atom.structure = structure[key]
#
#      # Some more post processing.
#      CA_list = [] # first build a list of the indices of all the alpha carbons
#      for atom_index in self.AllAtoms:
#        atom = self.AllAtoms[atom_index]
#        if (atom.residue.strip() in self.protein_resnames
#          and atom.atomname.strip() == "CA"):
#          CA_list.append(atom_index)
#
#      # some more post processing.
#      change = True
#      while change == True:
#        change = False
#
#        # A residue of index i is only going to be in an alpha helix
#        # its CA is within 6 A of the CA of the residue i + 3
#        for CA_atom_index in CA_list:
#          CA_atom = self.AllAtoms[CA_atom_index]
#          if CA_atom.structure == "ALPHA":
#            # so it's in an alpha helix
#            another_alpha_is_close = False
#            for other_CA_atom_index in CA_list:
#              # so now compare that CA to all the other CA's
#              other_CA_atom = self.AllAtoms[other_CA_atom_index]
#              if other_CA_atom.structure == "ALPHA": # so it's also in an alpha helix
#                if other_CA_atom.resid - 3 == CA_atom.resid or other_CA_atom.resid + 3 == CA_atom.resid:
#                  # so this CA atom is one of the ones the first atom
#                  # might hydrogen bond with
#                  if other_CA_atom.coordinates.dist_to(CA_atom.coordinates) < 6.0:
#                    # so these two CA atoms are close enough together
#                    # that their residues are probably hydrogen bonded
#                    another_alpha_is_close = True
#                    break
#            if another_alpha_is_close == False:
#              self.set_structure_of_residue(CA_atom.chain, CA_atom.resid, "OTHER")
#              change = True
#
#        # Alpha helices are only alpha helices if they span at least 4
#        # residues (to wrap around and hydrogen bond). I'm going to
#        # require them to span at least 5 residues, based on
#        # examination of many structures.
#        for index_in_list in range(len(CA_list)-5):
#
#          index_in_pdb1 = CA_list[index_in_list]
#          index_in_pdb2 = CA_list[index_in_list+1]
#          index_in_pdb3 = CA_list[index_in_list+2]
#          index_in_pdb4 = CA_list[index_in_list+3]
#          index_in_pdb5 = CA_list[index_in_list+4]
#          index_in_pdb6 = CA_list[index_in_list+5]
#
#          atom1 = self.AllAtoms[index_in_pdb1]
#          atom2 = self.AllAtoms[index_in_pdb2]
#          atom3 = self.AllAtoms[index_in_pdb3]
#          atom4 = self.AllAtoms[index_in_pdb4]
#          atom5 = self.AllAtoms[index_in_pdb5]
#          atom6 = self.AllAtoms[index_in_pdb6]
#
#          if (atom1.resid + 1 == atom2.resid
#            and atom2.resid + 1 == atom3.resid
#            and atom3.resid + 1 == atom4.resid
#            and atom4.resid + 1 == atom5.resid
#            and atom5.resid + 1 == atom6.resid): # so they are sequential
#            if atom1.structure != "ALPHA" and atom2.structure == "ALPHA" and atom3.structure != "ALPHA":
#              self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
#              change = True
#            if atom2.structure != "ALPHA" and atom3.structure == "ALPHA" and atom4.structure != "ALPHA":
#              self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
#              change = True
#            if atom3.structure != "ALPHA" and atom4.structure == "ALPHA" and atom5.structure != "ALPHA":
#              self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
#              change = True
#            if atom4.structure != "ALPHA" and atom5.structure == "ALPHA" and atom6.structure != "ALPHA":
#              self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
#              change = True
#
#            if (atom1.structure != "ALPHA"
#              and atom2.structure == "ALPHA"
#              and atom3.structure == "ALPHA"
#              and atom4.structure != "ALPHA"):
#              self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
#              self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
#              change = True
#            if atom2.structure != "ALPHA" and atom3.structure == "ALPHA" and atom4.structure == "ALPHA" and atom5.structure != "ALPHA":
#              self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
#              self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
#              change = True
#            if (atom3.structure != "ALPHA"
#              and atom4.structure == "ALPHA"
#              and atom5.structure == "ALPHA"
#              and atom6.structure != "ALPHA"):
#              self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
#              self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
#              change = True
#
#            if (atom1.structure != "ALPHA"
#              and atom2.structure == "ALPHA"
#              and atom3.structure == "ALPHA"
#              and atom4.structure == "ALPHA"
#              and atom5.structure != "ALPHA"):
#              self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
#              self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
#              self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
#              change = True
#            if (atom2.structure != "ALPHA"
#              and atom3.structure == "ALPHA"
#              and atom4.structure == "ALPHA"
#              and atom5.structure == "ALPHA"
#              and atom6.structure != "ALPHA"):
#              self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
#              self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
#              self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
#              change = True
#
#            if (atom1.structure != "ALPHA"
#              and atom2.structure == "ALPHA"
#              and atom3.structure == "ALPHA"
#              and atom4.structure == "ALPHA"
#              and atom5.structure == "ALPHA"
#              and atom6.structure != "ALPHA"):
#              self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
#              self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
#              self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
#              self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
#              change = True
#
#        # now go through each of the BETA CA atoms. A residue is only
#        # going to be called a beta sheet if CA atom is within 6.0 A
#        # of another CA beta, same chain, but index difference > 2.
#        for CA_atom_index in CA_list:
#            CA_atom = self.AllAtoms[CA_atom_index]
#            if CA_atom.structure == "BETA":
#              # so it's in a beta sheet
#              another_beta_is_close = False
#              for other_CA_atom_index in CA_list:
#                if other_CA_atom_index != CA_atom_index:
#                  # so not comparing an atom to itself
#                  other_CA_atom = self.AllAtoms[other_CA_atom_index]
#                  if other_CA_atom.structure == "BETA":
#                    # so you're comparing it only to other BETA-sheet atoms
#                    if other_CA_atom.chain == CA_atom.chain:
#                      # so require them to be on the same chain. needed to indecies can be fairly compared
#                      if math.fabs(other_CA_atom.resid - CA_atom.resid) > 2:
#                        # so the two residues are not simply adjacent to each other on the chain
#                        if CA_atom.coordinates.dist_to(other_CA_atom.coordinates) < 6.0:
#                          # so these to atoms are close to each other
#                          another_beta_is_close = True
#                          break
#              if another_beta_is_close == False:
#                self.set_structure_of_residue(CA_atom.chain, CA_atom.resid, "OTHER")
#                change = True
#
#        # Now some more post-processing needs to be done. Do this
#        # again to clear up mess that may have just been created
#        # (single residue beta strand, for example)
#        # Beta sheets are usually at least 3 residues long
#
#        for index_in_list in range(len(CA_list)-3):
#
#          index_in_pdb1 = CA_list[index_in_list]
#          index_in_pdb2 = CA_list[index_in_list+1]
#          index_in_pdb3 = CA_list[index_in_list+2]
#          index_in_pdb4 = CA_list[index_in_list+3]
#
#          atom1 = self.AllAtoms[index_in_pdb1]
#          atom2 = self.AllAtoms[index_in_pdb2]
#          atom3 = self.AllAtoms[index_in_pdb3]
#          atom4 = self.AllAtoms[index_in_pdb4]
#
#          if (atom1.resid + 1 == atom2.resid and atom2.resid + 1 ==
#            atom3.resid and atom3.resid + 1 == atom4.resid):
#            # so they are sequential
#
#            if (atom1.structure != "BETA"
#              and atom2.structure == "BETA"
#              and atom3.structure != "BETA"):
#              self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
#              change = True
#            if (atom2.structure != "BETA"
#              and atom3.structure == "BETA"
#              and atom4.structure != "BETA"):
#              self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
#              change = True
#            if (atom1.structure != "BETA"
#              and atom2.structure == "BETA"
#              and atom3.structure == "BETA"
#              and atom4.structure != "BETA"):
#              self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
#              self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
#              change = True
#
#    def set_structure_of_residue(self, chain, resid, structure):
#      for atom_index in self.AllAtoms:
#        atom = self.AllAtoms[atom_index]
#        if atom.chain == chain and atom.resid == resid:
#          atom.structure = structure

# TODO(bramsundar): Rework this to use numpy instead of custom
# implementations.
class MathFunctions:

#  def vector_subtraction(self, vector1, vector2): # vector1 - vector2
#    return point(vector1.x - vector2.x, vector1.y - vector2.y, vector1.z - vector2.z)
#
#  def CrossProduct(self, Pt1, Pt2): # never tested
#    Response = point(0,0,0)
#
#    Response.x = Pt1.y * Pt2.z - Pt1.z * Pt2.y
#    Response.y = Pt1.z * Pt2.x - Pt1.x * Pt2.z
#    Response.z = Pt1.x * Pt2.y - Pt1.y * Pt2.x
#
#    return Response;
#
#  def vector_scalar_multiply(self, vector, scalar):
#    return point(vector.x * scalar, vector.y * scalar, vector.z * scalar)
#
#  def dot_product(self, point1, point2):
#    return point1.x * point2.x + point1.y * point2.y + point1.z * point2.z
#
#  def dihedral(self, point1, point2, point3, point4): # never tested
#
#    b1 = self.vector_subtraction(point2, point1)
#    b2 = self.vector_subtraction(point3, point2)
#    b3 = self.vector_subtraction(point4, point3)
#
#    b2Xb3 = self.CrossProduct(b2,b3)
#    b1Xb2 = self.CrossProduct(b1,b2)
#
#    b1XMagb2 = self.vector_scalar_multiply(b1,b2.Magnitude())
#    radians = math.atan2(self.dot_product(b1XMagb2,b2Xb3), self.dot_product(b1Xb2,b2Xb3))
#    return radians
#
#  def angle_between_three_points(self, point1, point2, point3): # As in three connected atoms
#    vector1 = self.vector_subtraction(point1, point2)
#    vector2 = self.vector_subtraction(point3, point2)
#    return self.angle_between_points(vector1, vector2)
#
#  def angle_between_points(self, point1, point2):
#    new_point1 = self.return_normalized_vector(point1)
#    new_point2 = self.return_normalized_vector(point2)
#    dot_prod = self.dot_product(new_point1, new_point2)
#    if dot_prod > 1.0: dot_prod = 1.0 # to prevent errors that can rarely occur
#    if dot_prod < -1.0: dot_prod = -1.0
#    return math.acos(dot_prod)
#
#  def return_normalized_vector(self, vector):
#    dist = self.distance(point(0,0,0), vector)
#    return point(vector.x/dist, vector.y/dist, vector.z/dist)
#
  def distance(self, point1, point2):
    return point1.dist_to(point2)
#
#  def project_point_onto_plane(self, apoint, plane_coefficients):
#    # essentially finds the point on the plane that is closest to the
#    # specified point the plane_coefficients are [a,b,c,d], where the
#    # plane is ax + by + cz = d
#
#    # First, define a plane using cooeficients a, b, c, d such that ax + by + cz = d
#    a = plane_coefficients[0]
#    b = plane_coefficients[1]
#    c = plane_coefficients[2]
#    d = plane_coefficients[3]
#
#    # Now, define a point in space (s,u,v)
#    s = apoint.x
#    u = apoint.y
#    v = apoint.z
#
#    # the formula of a line perpendicular to the plan passing through (s,u,v) is:
#    #x = s + at
#    #y = u + bt
#    #z = v + ct
#
#    t = (d - a*s - b*u - c*v) / (a*a + b*b + c*c)
#
#    # here's the point closest on the plane
#    x = s + a*t
#    y = u + b*t
#    z = v + c*t
#
#    return point(x,y,z)
