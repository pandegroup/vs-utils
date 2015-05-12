"""
Helper Classes and Functions for docking fingerprint computation.

TODO(bramsundar): Most of the files below cannot be meaningfully tested
without a body of PDBs that exhibit the various amino-acids in questions.
Build up such a body of PDBs which we can use in the test work.

The code below contains heavily modified parts of Jacob Durrant's
NNScore 2.0.1. The following notice is copied from the original NNScore
file:
# NNScore 2.01 is released under the GNU General Public License (see
# http://www.gnu.org/licenses/gpl.html).
# If you have any questions, comments, or suggestions, please don't
# hesitate to contact me, Jacob Durrant, at jdurrant [at] ucsd [dot]
# edu. If you use NNScore 2.01 in your work, please cite [REFERENCE
# HERE].
"""
__author__ = "Bharath Ramsundar and Jacob Durrant"
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
                coordinates=Point(coords=np.array([99999, 99999, 99999])), element="",
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

  def side_chain_or_backbone(self):
    """Determine whether receptor atom belongs to residue sidechain or backbone.
    """
    # TODO(rbharath): Should this be an atom function?
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
