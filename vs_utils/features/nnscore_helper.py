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

  TODO(bramsundar): Ligand or Protein atom?
  """

  def __init__ (self):
    self.atomname = ""
    self.residue = ""
    self.coordinates = point(99999, 99999, 99999)
    self.element = ""
    self.PDBIndex = ""
    self.line=""
    self.atomtype=""
    self.IndicesOfAtomsConnecting=[]
    self.charge = 0
    self.resid = 0
    self.chain = ""
    self.structure = ""
    self.comment = ""

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
      #elif two_letters=='HG':
      #    self.element='HG'
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

## TODO(bramsundar): Uncomment this class once tests are written for
## classes above.
#class PDB:
#
#  def __init__ (self):
#    self.AllAtoms={}
#    self.NonProteinAtoms = {}
#    self.max_x = -9999.99
#    self.min_x = 9999.99
#    self.max_y = -9999.99
#    self.min_y = 9999.99
#    self.max_z = -9999.99
#    self.min_z = 9999.99
#    self.rotateable_bonds_count = 0
#    self.functions = MathFunctions()
#    self.protein_resnames = ["ALA", "ARG", "ASN", "ASP", "ASH", "ASX",
#      "CYS", "CYM", "CYX", "GLN", "GLU", "GLH", "GLX", "GLY", "HIS",
#      "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE",
#      "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
#    self.aromatic_rings = []
#    self.charges = [] # a list of points
#    self.OrigFileName = ""
#
#  def LoadPDB_from_file(self, FileName, line_header=""):
#
#    self.line_header=line_header
#
#    # Now load the file into a list
#    file = open(FileName,"r")
#    lines = file.readlines()
#    file.close()
#    self.LoadPDB_from_list(lines, self.line_header)
#
#  def LoadPDB_from_list(self, lines, line_header=""):
#
#    self.line_header=line_header
#    autoindex = 1
#    self.__init__()
#    # going to keep track of atomname_resid_chain pairs, to make sure
#    # redundants aren't loaded.  This basically gets rid of rotomers,
#    # I think.
#    atom_already_loaded = []
#
#    for t in range(0,len(lines)):
#      line=lines[t]
#
#      if "between atoms" in line and " A " in line:
#        self.rotateable_bonds_count = self.rotateable_bonds_count + 1
#
#      if len(line) >= 7:
#        # Load atom data (coordinates, etc.)
#        if line[0:4]=="ATOM" or line[0:6]=="HETATM":
#          TempAtom = atom()
#          TempAtom.ReadPDBLine(line)
#
#          # this string unique identifies each atom
#          key = (TempAtom.atomname.strip() + "_" +
#            str(TempAtom.resid) + "_" + TempAtom.residue.strip() +
#            "_" + TempAtom.chain.strip())
#          # so this is a receptor atom that has already been loaded once
#          if (key in atom_already_loaded and
#            TempAtom.residue.strip() in self.protein_resnames):
#            print (self.line_header + "WARNING: Duplicate receptor atom detected: \""
#             + TempAtom.line.strip()+ "\". Not loading this duplicate."
#
#          # so either the atom hasn't been loaded, or else it's a non-receptor
#          # atom so note that non-receptor atoms can have redundant names, but
#          # receptor atoms cannot.  This is because protein residues often
#          # contain rotamers
#          if (not key in atom_already_loaded or not
#            TempAtom.residue.strip() in self.protein_resnames):
#            # so each atom can only be loaded once. No rotamers.
#            atom_already_loaded.append(key)
#            # So you're actually reindexing everything here.
#            self.AllAtoms[autoindex] = TempAtom
#            if (not TempAtom.residue[-3:] in self.protein_resnames):
#              self.NonProteinAtoms[autoindex] = TempAtom
#
#            autoindex = autoindex + 1
#
#    self.CheckProteinFormat()
#    # only for the ligand, because bonds can be inferred based on
#    # atomnames from PDB
#    self.CreateBondsByDistance()
#      self.assign_aromatic_rings()
#      self.assign_charges()
#
#    def printout(self, thestring):
#      lines = textwrap.wrap(thestring, 80)
#      for line in lines:
#        print line
#
#    def SavePDB(self, filename):
#      f = open(filename, 'w')
#      towrite = self.SavePDBString()
#      # just so no PDB is empty, VMD will load them all
#      if towrite.strip() == "":
#        towrite = "ATOM      1  X   XXX             0.000   0.000   0.000                       X"
#      f.write(towrite)
#      f.close()
#
#    def SavePDBString(self):
#
#      ToOutput = ""
#
#      # write coordinates
#      for atomindex in self.AllAtoms:
#        ToOutput = ToOutput + self.AllAtoms[atomindex].CreatePDBLine(atomindex) + "\n"
#
#      return ToOutput
#
#    def AddNewAtom(self, atom):
#
#      # first get available index
#      t = 1
#      while t in self.AllAtoms.keys():
#        t = t + 1
#
#      # now add atom
#      self.AllAtoms[t] = atom
#
#    def connected_atoms_of_given_element(self, index, connected_atom_element):
#      atom = self.AllAtoms[index]
#      connected_atoms = []
#      for index2 in atom.IndeciesOfAtomsConnecting:
#        atom2 = self.AllAtoms[index2]
#        if atom2.element == connected_atom_element:
#          connected_atoms.append(index2)
#      return connected_atoms
#
#    def connected_heavy_atoms(self, index):
#      atom = self.AllAtoms[index]
#      connected_atoms = []
#      for index2 in atom.IndeciesOfAtomsConnecting:
#        atom2 = self.AllAtoms[index2]
#        if atom2.element != "H": connected_atoms.append(index2)
#      return connected_atoms
#
#    def CheckProteinFormat(self):
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
#          self.CheckProteinFormat_process_residue(residue, last_key)
#
#          residue = []
#          curr_res = key
#
#        residue.append(atom.atomname.strip())
#        last_key = key
#
#      self.CheckProteinFormat_process_residue(residue, last_key)
#
#
#    def CheckProteinFormat_process_residue(self, residue, last_key):
#      temp = last_key.strip().split("_")
#      resname = temp[0]
#      real_resname = resname[-3:]
#      resid = temp[1]
#      chain = temp[2]
#
#      if real_resname in self.protein_resnames: # so it's a protein residue
#
#        # TODO(bramsundar): The print statements here are highly
#        # redundant. Factor out the printing to avoid this block of
#        # strings.
#        if not "N" in residue:
#          self.printout('WARNING: There is no atom named "N" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine secondary structure. If this residue is far from the active site, this warning may not affect the NNScore.')
#          print ""
#        if not "C" in residue:
#          self.printout('WARNING: There is no atom named "C" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine secondary structure. If this residue is far from the active site, this warning may not affect the NNScore.')
#          print ""
#        if not "CA" in residue:
#          self.printout('WARNING: There is no atom named "CA" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine secondary structure. If this residue is far from the active site, this warning may not affect the NNScore.')
#          print ""
#
#        if real_resname == "GLU" or real_resname == "GLH" or real_resname == "GLX":
#          if not "OE1" in residue:
#            self.printout('WARNING: There is no atom named "OE1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "OE2" in residue:
#            self.printout('WARNING: There is no atom named "OE2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#
#        if real_resname == "ASP" or real_resname == "ASH" or real_resname == "ASX":
#          if not "OD1" in residue:
#            self.printout('WARNING: There is no atom named "OD1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#           not "OD2" in residue:
#            self.printout('WARNING: There is no atom named "OD2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#
#        if real_resname == "LYS" or real_resname == "LYN":
#          if not "NZ" in residue:
#            self.printout('WARNING: There is no atom named "NZ" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#
#        if real_resname == "ARG":
#          if not "NH1" in residue:
#            self.printout('WARNING: There is no atom named "NH1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "NH2" in residue:
#            self.printout('WARNING: There is no atom named "NH2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#
#        if real_resname == "HIS" or real_resname == "HID" or real_resname == "HIE" or real_resname == "HIP":
#          if not "NE2" in residue:
#            self.printout('WARNING: There is no atom named "NE2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "ND1" in residue:
#            self.printout('WARNING: There is no atom named "ND1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#
#        if real_resname == "PHE":
#          if not "CG" in residue:
#            self.printout('WARNING: There is no atom named "CG" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CD1" in residue:
#            self.printout('WARNING: There is no atom named "CD1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CD2" in residue:
#            self.printout('WARNING: There is no atom named "CD2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CE1" in residue:
#            self.printout('WARNING: There is no atom named "CE1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CE2" in residue:
#            self.printout('WARNING: There is no atom named "CE2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CZ" in residue:
#            self.printout('WARNING: There is no atom named "CZ" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#
#        if real_resname == "TYR":
#          if not "CG" in residue:
#            self.printout('WARNING: There is no atom named "CG" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CD1" in residue:
#            self.printout('WARNING: There is no atom named "CD1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CD2" in residue:
#            self.printout('WARNING: There is no atom named "CD2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CE1" in residue:
#            self.printout('WARNING: There is no atom named "CE1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CE2" in residue:
#            self.printout('WARNING: There is no atom named "CE2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CZ" in residue:
#            self.printout('WARNING: There is no atom named "CZ" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#
#        if real_resname == "TRP":
#          if not "CG" in residue:
#            self.printout('WARNING: There is no atom named "CG" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CD1" in residue:
#            self.printout('WARNING: There is no atom named "CD1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CD2" in residue:
#            self.printout('WARNING: There is no atom named "CD2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "NE1" in residue:
#            self.printout('WARNING: There is no atom named "NE1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CE2" in residue:
#            self.printout('WARNING: There is no atom named "CE2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CE3" in residue:
#            self.printout('WARNING: There is no atom named "CE3" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CZ2" in residue:
#            self.printout('WARNING: There is no atom named "CZ2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CZ3" in residue:
#            self.printout('WARNING: There is no atom named "CZ3" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CH2" in residue:
#            self.printout('WARNING: There is no atom named "CH2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#
#        if (real_resname == "HIS" or real_resname == "HID" or
#          real_resname == "HIE" or real_resname == "HIP"):
#          if not "CG" in residue:
#            self.printout('WARNING: There is no atom named "CG" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "ND1" in residue:
#            self.printout('WARNING: There is no atom named "ND1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CD2" in residue:
#            self.printout('WARNING: There is no atom named "CD2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "CE1" in residue:
#            self.printout('WARNING: There is no atom named "CE1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#          if not "NE2" in residue:
#            self.printout('WARNING: There is no atom named "NE2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
#            print ""
#
#
#    # Functions to determine the bond connectivity based on distance
#    # ==============================================================
#
#    def CreateBondsByDistance(self):
#      for AtomIndex1 in self.NonProteinAtoms:
#        atom1 = self.NonProteinAtoms[AtomIndex1]
#        if not atom1.residue[-3:] in self.protein_resnames: # so it's not a protein residue
#          for AtomIndex2 in self.NonProteinAtoms:
#            if AtomIndex1 != AtomIndex2:
#              atom2 = self.NonProteinAtoms[AtomIndex2]
#              if not atom2.residue[-3:] in self.protein_resnames: # so it's not a protein residue
#                dist = self.functions.distance(atom1.coordinates, atom2.coordinates)
#
#                if (dist < self.BondLength(atom1.element, atom2.element) * 1.2):
#                  atom1.AddNeighborAtomIndex(AtomIndex2)
#                  atom2.AddNeighborAtomIndex(AtomIndex1)
#
#    def BondLength(self, element1, element2):
#      """Bond lengths taken from Handbook of Chemistry and Physics.
#      The information provided there was very specific, so I tried
#      to pick representative examples and used the bond lengths from
#      those. Sitautions could arise where these lengths would be
#      incorrect, probably slight errors (<0.06) in the hundreds."""
#
#      distance = 0.0
#      if element1 == "C" and element2 == "C":
#        distance = 1.53
#      if element1 == "N" and element2 == "N":
#        distance = 1.425
#      if element1 == "O" and element2 == "O":
#        distance = 1.469
#      if element1 == "S" and element2 == "S":
#        distance = 2.048
#      if ((element1 == "C" and element2 == "H")
#        or (element1 == "H" and element2 == "C")):
#        distance = 1.059
#      if ((element1 == "C" and element2 == "N")
#        or (element1 == "N" and element2 == "C")):
#        distance = 1.469
#      if ((element1 == "C" and element2 == "O")
#        or (element1 == "O" and element2 == "C")):
#        distance = 1.413
#      if ((element1 == "C" and element2 == "S")
#        or (element1 == "S" and element2 == "C")):
#        distance = 1.819
#      if ((element1 == "N" and element2 == "H")
#        or (element1 == "H" and element2 == "N")):
#        distance = 1.009
#      if ((element1 == "N" and element2 == "O")
#        or (element1 == "O" and element2 == "N")):
#        distance = 1.463
#      if ((element1 == "O" and element2 == "S")
#        or (element1 == "S" and element2 == "O")):
#        distance = 1.577
#      if ((element1 == "O" and element2 == "H")
#        or (element1 == "H" and element2 == "O")):
#        distance = 0.967
#      if ((element1 == "S" and element2 == "H")
#        or (element1 == "H" and element2 == "S")):
#        # This one not from source sited above. Not sure where it's from, but
#        # it wouldn't ever be used in the current context ("AutoGrow")
#        distance = 2.025/1.5
#      if ((element1 == "S" and element2 == "N")
#        or (element1 == "N" and element2 == "S")):
#        distance = 1.633
#
#      if ((element1 == "C" and element2 == "F")
#        or (element1 == "F" and element2 == "C")):
#        distance = 1.399
#      if ((element1 == "C" and element2 == "CL")
#        or (element1 == "CL" and element2 == "C")):
#        distance = 1.790
#      if ((element1 == "C" and element2 == "BR")
#        or (element1 == "BR" and element2 == "C")):
#        distance = 1.910
#      if ((element1 == "C" and element2 == "I")
#        or (element1 == "I" and element2 == "C")):
#        distance=2.162
#
#      if ((element1 == "S" and element2 == "BR")
#        or (element1 == "BR" and element2 == "S")):
#        distance = 2.321
#      if ((element1 == "S" and element2 == "CL")
#        or (element1 == "CL" and element2 == "S")):
#        distance = 2.283
#      if ((element1 == "S" and element2 == "F")
#        or (element1 == "F" and element2 == "S")):
#        distance = 1.640
#      if ((element1 == "S" and element2 == "I")
#        or (element1 == "I" and element2 == "S")):
#        distance= 2.687
#
#      if ((element1 == "P" and element2 == "BR")
#        or (element1 == "BR" and element2 == "P")):
#        distance = 2.366
#      if ((element1 == "P" and element2 == "CL")
#        or (element1 == "CL" and element2 == "P")):
#        distance = 2.008
#      if ((element1 == "P" and element2 == "F")
#        or (element1 == "F" and element2 == "P")):
#        distance = 1.495
#      if ((element1 == "P" and element2 == "I")
#        or (element1 == "I" and element2 == "P")):
#        distance= 2.490
#      if ((element1 == "P" and element2 == "O")
#        or (element1 == "O" and element2 == "P")):
#        distance= 1.6 # estimate based on eye balling Handbook of Chemistry and Physics
#
#      if ((element1 == "N" and element2 == "BR")
#        or (element1 == "BR" and element2 == "N")):
#        distance = 1.843
#      if ((element1 == "N" and element2 == "CL")
#        or (element1 == "CL" and element2 == "N")):
#        distance = 1.743
#      if ((element1 == "N" and element2 == "F")
#        or (element1 == "F" and element2 == "N")):
#        distance = 1.406
#      if ((element1 == "N" and element2 == "I")
#        or (element1 == "I" and element2 == "N")):
#        distance= 2.2
#
#      if ((element1 == "SI" and element2 == "BR")
#        or (element1 == "BR" and element2 == "SI")):
#        distance = 2.284
#      if ((element1 == "SI" and element2 == "CL")
#        or (element1 == "CL" and element2 == "SI")):
#        distance = 2.072
#      if ((element1 == "SI" and element2 == "F")
#        or (element1 == "F" and element2 == "SI")):
#        distance = 1.636
#      if ((element1 == "SI" and element2 == "P")
#        or (element1 == "P" and element2 == "SI")):
#        distance= 2.264
#      if ((element1 == "SI" and element2 == "S")
#        or (element1 == "S" and element2 == "SI")):
#        distance= 2.145
#      if ((element1 == "SI" and element2 == "SI")
#        or (element1 == "SI" and element2 == "SI")):
#        distance= 2.359
#      if ((element1 == "SI" and element2 == "C")
#        or (element1 == "C" and element2 == "SI")):
#        distance= 1.888
#      if ((element1 == "SI" and element2 == "N")
#        or (element1 == "N" and element2 == "SI")):
#        distance= 1.743
#      if ((element1 == "SI" and element2 == "O")
#        or (element1 == "O" and element2 == "SI")):
#        distance= 1.631
#
#      return distance;
#
#    # Functions to identify positive charges
#    # ======================================
#
#    def assign_charges(self):
#      # Get all the quartinary amines on non-protein residues (these
#      # are the only non-protein groups that will be identified as
#      # positively charged)
#      AllCharged = []
#      for atom_index in self.NonProteinAtoms:
#          atom = self.NonProteinAtoms[atom_index]
#          if (atom.element == "MG" or atom.element == "MN" or
#          atom.element == "RH" or atom.element == "ZN" or
#          atom.element == "FE" or atom.element == "BI" or
#          atom.element == "AS" or atom.element == "AG"):
#            chrg = self.charged(atom.coordinates, [atom_index], True)
#            self.charges.append(chrg)
#
#          if atom.element == "N":
#            # a quartinary amine, so it's easy
#            if atom.NumberOfNeighbors() == 4:
#              indexes = [atom_index]
#              indexes.extend(atom.IndeciesOfAtomsConnecting)
#              # so the indicies stored is just the index of the nitrogen and any attached atoms
#              chrg = self.charged(atom.coordinates, indexes, True)
#              self.charges.append(chrg)
#            # maybe you only have two hydrogens added, by they're sp3 hybridized.
#            # Just count this as a quartinary amine, since I think the positive
#            # charge would be stabalized.
#            elif atom.NumberOfNeighbors() == 3:
#              nitrogen = atom
#              atom1 = self.AllAtoms[atom.IndeciesOfAtomsConnecting[0]]
#              atom2 = self.AllAtoms[atom.IndeciesOfAtomsConnecting[1]]
#              atom3 = self.AllAtoms[atom.IndeciesOfAtomsConnecting[2]]
#              angle1 = (self.functions.angle_between_three_points(atom1.coordinates,
#                nitrogen.coordinates, atom2.coordinates) * 180.0 /
#                math.pi)
#              angle2 = (self.functions.angle_between_three_points(atom1.coordinates,
#                nitrogen.coordinates, atom3.coordinates) * 180.0 /
#                math.pi)
#              angle3 = (self.functions.angle_between_three_points(atom2.coordinates,
#                nitrogen.coordinates, atom3.coordinates) * 180.0 /
#                math.pi)
#              average_angle = (angle1 + angle2 + angle3) / 3
#              if math.fabs(average_angle - 109.0) < 5.0:
#                indexes = [atom_index]
#                indexes.extend(atom.IndeciesOfAtomsConnecting)
#                # so indexes added are the nitrogen and any attached atoms.
#                chrg = self.charged(nitrogen.coordinates, indexes, True)
#                self.charges.append(chrg)
#
#          # let's check for guanidino-like groups (actually H2N-C-NH2,
#          # where not CN3.)
#          if atom.element == "C":
#            # the carbon has only three atoms connected to it
#            if atom.NumberOfNeighbors() == 3:
#              nitrogens = self.connected_atoms_of_given_element(atom_index,"N")
#              # so carbon is connected to at least two nitrogens now
#              # we need to count the number of nitrogens that are only
#              # connected to one heavy atom (the carbon)
#              if len(nitrogens) >= 2:
#
#                nitrogens_to_use = []
#                all_connected = atom.IndeciesOfAtomsConnecting[:]
#                not_isolated = -1
#
#                for atmindex in nitrogens:
#                  if len(self.connected_heavy_atoms(atmindex)) == 1:
#                      nitrogens_to_use.append(atmindex)
#                      all_connected.remove(atmindex)
#
#                if len(all_connected) > 0:
#                  # get the atom that connects this charged group to
#                  # the rest of the molecule, ultimately to make sure
#                  # it's sp3 hybridized
#                  not_isolated = all_connected[0]
#
#                # so there are at two nitrogens that are only
#                # connected to the carbon (and probably some
#                # hydrogens)
#                if len(nitrogens_to_use) == 2 and not_isolated != -1:
#
#                  # now you need to make sure not_isolated atom is sp3 hybridized
#                  not_isolated_atom = self.AllAtoms[not_isolated]
#                  if ((not_isolated_atom.element == "C" and
#                      not_isolated_atom.NumberOfNeighbors()==4)
#                    or (not_isolated_atom.element == "O"
#                      and not_isolated_atom.NumberOfNeighbors()==2)
#                    or not_isolated_atom.element == "N"
#                    or not_isolated_atom.element == "S"
#                    or not_isolated_atom.element == "P"):
#
#                    pt = self.AllAtoms[nitrogens_to_use[0]].coordinates.copy_of()
#                    pt.x = pt.x + self.AllAtoms[nitrogens_to_use[1]].coordinates.x
#                    pt.y = pt.y + self.AllAtoms[nitrogens_to_use[1]].coordinates.y
#                    pt.z = pt.z + self.AllAtoms[nitrogens_to_use[1]].coordinates.z
#                    pt.x = pt.x / 2.0
#                    pt.y = pt.y / 2.0
#                    pt.z = pt.z / 2.0
#
#                    indexes = [atom_index]
#                    indexes.extend(nitrogens_to_use)
#                    indexes.extend(self.connected_atoms_of_given_element(nitrogens_to_use[0],"H"))
#                    indexes.extend(self.connected_atoms_of_given_element(nitrogens_to_use[1],"H"))
#
#                    chrg = self.charged(pt, indexes, True) # True because it's positive
#                    self.charges.append(chrg)
#
#          if atom.element == "C": # let's check for a carboxylate
#              # a carboxylate carbon will have three items connected to it.
#              if atom.NumberOfNeighbors() == 3:
#                oxygens = self.connected_atoms_of_given_element(atom_index,"O")
#                # a carboxylate will have two oxygens connected to
#                # it. Now, each of the oxygens should be connected
#                # to only one heavy atom (so if it's connected to a
#                # hydrogen, that's okay)
#                if len(oxygens) == 2:
#                  if (len(self.connected_heavy_atoms(oxygens[0])) == 1
#                    and len(self.connected_heavy_atoms(oxygens[1])) == 1):
#                    # so it's a carboxylate! Add a negative charge.
#                    pt = self.AllAtoms[oxygens[0]].coordinates.copy_of()
#                    pt.x = pt.x + self.AllAtoms[oxygens[1]].coordinates.x
#                    pt.y = pt.y + self.AllAtoms[oxygens[1]].coordinates.y
#                    pt.z = pt.z + self.AllAtoms[oxygens[1]].coordinates.z
#                    pt.x = pt.x / 2.0
#                    pt.y = pt.y / 2.0
#                    pt.z = pt.z / 2.0
#                    chrg = self.charged(pt, [oxygens[0],
#                        atom_index, oxygens[1]], False)
#                    self.charges.append(chrg)
#
#          # let's check for a phosphate or anything where a phosphorus is bound
#          # to two oxygens where both oxygens are bound to only one heavy atom
#          # (the phosphorus). I think this will get several phosphorus
#          # substances.
#          if atom.element == "P":
#            oxygens = self.connected_atoms_of_given_element(atom_index,"O")
#            if len(oxygens) >=2: # the phosphorus is bound to at least two oxygens
#              # now count the number of oxygens that are only bound to the phosphorus
#              count = 0
#              for oxygen_index in oxygens:
#                if len(self.connected_heavy_atoms(oxygen_index)) == 1: count = count + 1
#              if count >=2: # so there are at least two oxygens that are only bound to the phosphorus
#                indexes = [atom_index]
#                indexes.extend(oxygens)
#                chrg = self.charged(atom.coordinates, indexes, False)
#                self.charges.append(chrg)
#
#          # let's check for a sulfonate or anything where a sulfur is
#          # bound to at least three oxygens and at least three are
#          # bound to only the sulfur (or the sulfur and a hydrogen).
#          if atom.element == "S":
#            oxygens = self.connected_atoms_of_given_element(atom_index,"O")
#            # the sulfur is bound to at least three oxygens now
#            # count the number of oxygens that are only bound to the
#            # sulfur
#            if len(oxygens) >=3:
#              count = 0
#              for oxygen_index in oxygens:
#                if len(self.connected_heavy_atoms(oxygen_index)) == 1: count = count + 1
#              # so there are at least three oxygens that are only
#              # bound to the sulfur
#              if count >=3:
#                indexes = [atom_index]
#                indexes.extend(oxygens)
#                chrg = self.charged(atom.coordinates, indexes, False)
#                self.charges.append(chrg)
#
#      # Now that you've found all the positive charges in non-protein
#      # residues, it's time to look for aromatic rings in protein
#      # residues
#      curr_res = ""
#      first = True
#      residue = []
#
#      for atom_index in self.AllAtoms:
#        atom = self.AllAtoms[atom_index]
#        key = atom.residue + "_" + str(atom.resid) + "_" + atom.chain
#
#        if first == True:
#          curr_res = key
#          first = False
#
#        if key != curr_res:
#
#          self.assign_charged_from_protein_process_residue(residue, last_key)
#          residue = []
#          curr_res = key
#
#        residue.append(atom_index)
#        last_key = key
#
#      self.assign_charged_from_protein_process_residue(residue, last_key)
#
#    def assign_charged_from_protein_process_residue(self, residue, last_key):
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
#            chrg = self.charged(atom.coordinates, indexes, True)
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
#            chrg = self.charged(charge_pt, indices, True)
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
#              chrg = self.charged(charge_pt, indices, True)
#              self.charges.append(chrg)
#
#        if real_resname == "GLU" or real_resname == "GLH" or real_resname == "GLX":
#          # regardless of protonation state, assume it's charged. This based on
#          # "The Cation-Pi Interaction," which suggests protonated state would
#          # be stabalized.
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
#            chrg = self.charged(charge_pt, indices, False) # False because it's a negative charge
#            self.charges.append(chrg)
#
#        # TODO(bramsundar): This comment about Cation-Pi interactions
#        # is repeated in multiple places. Look into this interaction
#        # and verify that it holds true for the residues in question.
#        if (real_resname == "ASP" or real_resname == "ASH" or
#          real_resname == "ASX"):
#          # regardless of protonation state, assume it's charged. This based on
#          # "The Cation-Pi Interaction," which suggests protonated state would
#          # be stabalized.
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
#              chrg = self.charged(charge_pt, indices, False)
#              self.charges.append(chrg)
#
#    class charged():
#      def __init__(self, coordinates, indices, positive):
#        self.coordinates = coordinates
#        self.indices = indices
#        # true or false to specifiy if positive or negative charge
#        self.positive = positive
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
#          for substituent_atom_index in cur_atom.IndeciesOfAtomsConnecting:
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
#      for conneceted_atom in atom.IndeciesOfAtomsConnecting:
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
#      for conneceted_atom in atom.IndeciesOfAtomsConnecting:
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
