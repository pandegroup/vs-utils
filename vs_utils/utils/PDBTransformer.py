#Written by Evan N. Feinberg at Stanford University, contact: enf@stanford.edu

'''
This code is temporarily deprecated and will ultimately be removed. If it's still here by 12/1/2015, 
please yell at enf@stanford.edu
'''


from nnscore_pdb import PDB 
from copy import deepcopy
import numpy as np
import pickle 

'''Pseudocode 

1. Create objects of class PDB for protein and ligand
2. Compute centroid of ligand
3. translate coordinates of all other atoms such that ligand centroid is now origin, for both 
  ligand PDB and protein PDB
4. combine PDB objects
5. Write out new PDB (.pdb and pickle of PDB() object)
'''

'''
Pseudocode, detail: 

1. Create objects of class PDB for protein and ligand 
2. Compute centroid of ligand
  a. define a function for computing centroid given a PDB object. Maybe extend to make it 
    accept just coordinates?
    -further decomposition:
      Take object as input --> convert to coordinates --> then send to separate method 
      that computes centroid 
3. For each atom for each PDB oject: 
  -change x, y, and z coordinates by subtracting out centroid x, y, and z
4. Function that combines PDB objects: 
  def mergePDBs(pdb1, pdb2)
    copy pdb1
    add all from pdb2 to pdb1, but offset atom indices by max(atom index of pdb 1)
'''


def compute_centroid(molecule):
  '''given molecule, an instance of class PDB, compute the x,y,z centroid of that molecule'''

  coordinates = [molecule.all_atoms[atom].coordinates.coords for atom in molecule.all_atoms]
  centroid = np.sum(coordinates, axis=0)/float(len(coordinates))
  return(centroid)

def generate_random_unit_vector():
  '''generate a random unit vector on the 3-sphere
  citation:
  http://mathworld.wolfram.com/SpherePointPicking.html

  a. Choose random theta \element [0, 2*pi]
  b. Choose random z \element [-1, 1]
  c. Compute output: (x,y,z) = (sqrt(1-z^2)*cos(theta), sqrt(1-z^2)*sin(theta),z)
  d. output u

  '''

  theta = np.random.uniform(low = 0.0, high = 2*np.pi)
  z = np.random.uniform(low = -1.0, high = 1.0)
  u = np.array([np.sqrt(1-z**2)*np.cos(theta), np.sqrt(1-z**2)*np.sin(theta), z])
  return(u)

def generate_random_rotation_matrix():
  '''
   1. generate a random unit vector, i.e., randomly sampled from the unit 3-sphere
      a. see function generate_random_unit_vector() for details
    2. Generate a second random unit vector thru the algorithm in (1), output v
      a. If absolute value of u \dot v > 0.99, repeat. This is important for numerical stability 
        (intuition: we want them to be as linearly independent as possible or else the 
          orthogonalized version of v will be much shorter in magnitude compared to u. I assume
          in Stack they took this from Gram-Schmidt orthogonalization?)
      b. v' = v - (u \dot v)*u, i.e. subtract out the component of v that's in u's direction 
      c. normalize v' (this isn't in Stack but I assume it must be done)
    3. find w = u \cross v' 
    4. u, v', and w will form the columns of a rotation matrix, R. The intuition is that
        u, v' and w are, respectively, what the standard basis vectors e1, e2, and e3 will be 
        mapped to under the transformation.
  '''

  u = generate_random_unit_vector()
  v = generate_random_unit_vector()
  while np.abs(np.dot(u,v)) >= 0.99:
    v = generate_random_unit_vector()

  vp = v - (np.dot(u,v)*u)
  vp /= np.linalg.norm(vp)

  w = np.cross(u,vp)

  R = np.column_stack((u, vp, w))
  return(R)


class PDBTransformer:
  def __init__(self):
    self.protein = None
    self.ligand = None
    self.system = None
    self.box = None
    self.box_x = 100000.0
    self.box_y = 100000.0
    self.box_z = 100000.0


  def transform(self, protein_pdb, protein_pdbqt, ligand_pdb, ligand_pdbqt, 
                save_dir, box_x, box_y, box_z, nb_rotations = 0, nb_reflections=0):
    '''Takes as input files (strings) for pdb/pdbqt of the protein, pdb/pdbqt of the ligand, a filename to be saved 
    of the merged system called system_pdb, a filename to be saved of the box called box_pdb, a filename with a .p extension
    of the box in pickle format called box_pickle, and 3 floats for the x,y,z dimensions of the box

    This function then computes the centroid of the ligand; decrements this centroid from the atomic coordinates of protein and
    ligand atoms, and then merges the translated protein and ligand. This combined system/complex is then saved.

    This function then removes all atoms outside of the given box dimensions, and saves the resulting box in PDB file format
    as well as in a pickle format of the box's instance of the PDB object.
    '''


    self.protein = PDB()
    self.protein.load_from_files(protein_pdb, protein_pdbqt)

    self.ligand = PDB()
    self.ligand.load_from_files(ligand_pdb, ligand_pdbqt)

    self.box_x = float(box_x)
    self.box_y = float(box_y)
    self.box_z = float(box_z)

    protein_name = str(protein_pdb).split("/")[len(str(protein_pdb).split("/"))-2]
    system_pdb_file = "%s/%s.pdb" %(save_dir, protein_name)
    system_pickle_file = "%s/%s.pickle" %(save_dir, protein_name)

    ligand_centroid = compute_centroid(self.ligand)
    self.ligand = self.subtract_centroid(self.ligand, ligand_centroid)
    self.protein = self.subtract_centroid(self.protein, ligand_centroid)

    self.system = self.merge_molecules(self.protein, self.ligand)
    original_system = deepcopy(self.system)
    self.box = self.generate_box(original_system)
    original_box = deepcopy(self.box)

    transformed_systems = {}
    transformed_boxes = {}
    transformed_systems[(0,0)] = original_system
    transformed_boxes[(0,0)] = original_box

    for i in range(0,int(nb_rotations)):
      rotated_system = self.rotate_molecule(original_system)
      transformed_systems[(i+1,0)] =  rotated_system
      for j in range(0,int(nb_reflections)):
        reflected_system = self.reflect_molecule(rotated_system)
        transformed_systems[(i+1,j+1)] = reflected_system

    transformed_boxes = {}
    for key, transformed_system in transformed_systems.iteritems(): 
      transformed_box = self.generate_box(transformed_system)
      transformed_boxes[key] = transformed_box

    for key, transformed_system in transformed_systems.iteritems():
      print(key)
      print("%s/%s_%d_%d_system.pdb" %(save_dir, protein_name, key[0], key[1]))
      transformed_system.save_pdb("%s/%s_%d_%d_system.pdb" %(save_dir, protein_name, key[0], key[1]))
      pickle.dump(transformed_system, open("%s/%s_%d_%d_system.pickle" %(save_dir, protein_name, key[0], key[1]), "wb"))

    for key, transformed_box in transformed_boxes.iteritems():
      transformed_box.save_pdb("%s/%s_%d_%d_box.pdb" %(save_dir, protein_name, key[0], key[1]))
      pickle.dump(transformed_box, open("%s/%s_%d_%d_box.pickle" %(save_dir, protein_name, key[0], key[1]), "wb"))

    return(transformed_boxes)



  def generate_box(self, mol):
    '''
    generate_box takes as input a molecule of class PDB and removes all atoms outside of the given box dims
    '''
    molecule = deepcopy(mol)
    atoms_to_remove = []
    for index, atom in molecule.all_atoms.iteritems():
      coords = np.abs(atom.coordinates.coords)
      if coords[0] >= (self.box_x/2.) or coords[1] >= (self.box_y/2.) or coords[2] >= (self.box_z/2.):
        atoms_to_remove.append((index, atom))

    for index, atom in atoms_to_remove:
      molecule = self.remove_atom(molecule, index, atom)

    return(molecule)

  def remove_atom(self, molecule, atom_index, atom):
    '''
    remove_atom works by simply deleting the entry corresponding to that atom from the 
    molecule.all_atoms dictionary.
    '''


    if molecule.all_atoms[atom_index] != atom:
      print("Removing atoms from dictionary not safe!!")
      return

    del molecule.all_atoms[atom_index]
    return(molecule)

  def subtract_centroid(self, molecule, centroid):
    '''
    subtracts the centroid, a numpy array of dim 3, from all coordinates of all atoms in the molecule
    '''


    for atom in molecule.all_atoms:
      molecule.all_atoms[atom].coordinates.coords -= centroid
    return(molecule)

  def rotate_molecule(self, mol):
    '''
    Pseudocode:
    1. Generate random rotation matrix. This matrix applies a random transformation to any 
      3-vector such that, were the random transformation repeatedly applied, it would randomly
      sample along the surface of a sphere with radius equal to the norm of the given 3-vector
      cf. generate_random_rotation_matrix() for details
    2. Apply R to all atomic coordinatse. 
    3. Return rotated molecule
    '''
    molecule = deepcopy(mol)
    R = generate_random_rotation_matrix()
    all_coordinates = np.column_stack([molecule.all_atoms[index].coordinates.coords for index in molecule.all_atoms.iterkeys()])
    rotated_coordinates = np.dot(R, all_coordinates)
    for j, index in enumerate(molecule.all_atoms.iterkeys()):
      molecule.all_atoms[index].coordinates.coords = rotated_coordinates[:,j]

    return(molecule)

  def reflect_molecule(self, mol):
    '''
    Pseudocode:
    1. Generate unit vector that is randomly distributed around 3-sphere
    2. For each atom, project its coordinates onto the random unit vector from (1),
      and subtract twice the projection from the original coordinates to obtain its reflection 
    '''
    molecule = deepcopy(mol)
    a = generate_random_unit_vector()
    for index, atom in molecule.all_atoms.iteritems():
      v = atom.coordinates.coords 
      reflected_coords = v - 2. * (np.dot(a,v) / (np.dot(a,a)) * a)
      molecule.all_atoms[index].coordinates.coords = reflected_coords
    return(molecule)

  def merge_molecules(self, protein, ligand):
    '''
    Takes as input protein and ligand objects of class PDB and adds ligand atoms to the protein,
    and returns the new instance of class PDB called system that contains both sets of atoms.
    '''

    system = deepcopy(protein)
    greatest_index = len(protein.all_atoms) + 1
    autoindex = greatest_index + 1
    for index, ligand_atom in ligand.all_atoms.iteritems():
      system.add_new_atom(ligand_atom)
    return(system)


