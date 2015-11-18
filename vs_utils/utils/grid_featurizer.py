#Written by Evan N. Feinberg at Stanford University, contact: enf@stanford.edu
from copy import deepcopy
import numpy as np
import pickle
import gzip
import mdtraj as md
import time
from pybel import *
import networkx as nx
from collections import deque
import hashlib


def compute_centroid(coordinates):
  '''given molecule, an instance of class PDB, compute the x,y,z centroid of that molecule'''

  print(coordinates)
  centroid = np.mean(coordinates, axis=0)
  print(centroid)
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

def compute_pairwise_distances(protein_xyz, ligand_xyz):
  '''
  Takes an input m x 3 and n x 3 np arrays of 3d coords of protein and ligand,
  respectively, and outputs an m x n np array of pairwise distances in Angstroms
  between protein and ligand atoms. entry (i,j) is dist between the i'th protein
  atom and the j'th ligand atom
  '''

  pairwise_distances = np.zeros((np.shape(protein_xyz)[0], np.shape(ligand_xyz)[0]))
  for j in range(0, np.shape(ligand_xyz)[0]):
    differences = protein_xyz - ligand_xyz[j,:]
    squared_differences = np.square(differences)
    pairwise_distances[:,j] = np.sqrt(np.sum(squared_differences,1))

  return(pairwise_distances)

'''following two functions adapted from:
http://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
'''

def unit_vector(vector):
  """ Returns the unit vector of the vector.  """
  return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
  """ Returns the angle in radians between vectors 'v1' and 'v2'::

          >>> angle_between((1, 0, 0), (0, 1, 0))
          1.5707963267948966
          >>> angle_between((1, 0, 0), (1, 0, 0))
          0.0
          >>> angle_between((1, 0, 0), (-1, 0, 0))
          3.141592653589793
  """
  v1_u = unit_vector(v1)
  v2_u = unit_vector(v2)
  angle = np.arccos(np.dot(v1_u, v2_u))
  if np.isnan(angle):
    if (v1_u == v2_u).all():
      return 0.0
    else:
      return np.pi
  return angle

class grid_featurizer:
  def __init__(self, box_x=16.0, box_y=16.0, box_z=16.0, 
    nb_rotations = 0, nb_reflections=0, feature_types="ecfp", 
    ecfp_degree=2, ecfp_power = 8, splif_power = 8, 
    save_intermediates=False, ligand_only = False, **kwargs):

    self.box = None

    self.box_x = 100000.0
    self.box_y = 100000.0
    self.box_z = 100000.0

    self.box_x = float(box_x)/10.0
    self.box_y = float(box_y)/10.0
    self.box_z = float(box_z)/10.0


    self.ecfp_degree = ecfp_degree
    self.ecfp_power = ecfp_power
    self.splif_power = splif_power

    self.nb_rotations = nb_rotations
    self.nb_reflections = nb_reflections
    self.feature_types = feature_types

    self.save_intermediates = save_intermediates
    self.ligand_only = ligand_only

    self.hbond_angle_cutoff = 90.
    self.hbond_dist_cutoff = 4.0

  def transform(self, protein_pdb, ligand_pdb, save_dir):
    '''Takes as input files (strings) for pdb of the protein, pdb of the ligand, and a directory 
    to save intermediate files. 

    This function then computes the centroid of the ligand; decrements this centroid from the atomic coordinates of protein and
    ligand atoms, and then merges the translated protein and ligand. This combined system/complex is then saved.

    This function then computes a featurization with scheme specified by the user.
    '''
    a = time.time()
    protein_name = str(protein_pdb).split("/")[len(str(protein_pdb).split("/"))-2]
    system_pdb_file = "%s/%s.pdb" %(save_dir, protein_name)
    system_pickle_file = "%s/%s.pickle" %(save_dir, protein_name)

    if not self.ligand_only:
      protein_xyz, protein_ob = self.load_molecule(protein_pdb) 
    ligand_xyz, ligand_ob = self.load_molecule(ligand_pdb)

    if "ecfp" in self.feature_types:
      ecfp_array = self.compute_ecfp_features(ligand_ob, self.ecfp_degree, self.ecfp_power)
      with gzip.open("%s/%s_ecfp_array.pkl.gz" %(save_dir, protein_name), "wb") as f:
        pickle.dump(ecfp_array, f)
      return({(0,0): ecfp_array})

    centroid = compute_centroid(ligand_xyz)
    ligand_xyz = self.subtract_centroid(ligand_xyz, centroid)
    if not self.ligand_only:
      protein_xyz = self.subtract_centroid(protein_xyz, centroid)
    print(time.time()-a)

    if "splif" in self.feature_types:
      splif_begin = time.time()
      splif_array = self.featurize_splif(protein_xyz, protein_ob, ligand_xyz, ligand_ob)
      splif_end = time.time()
      print("Took %f seconds to featurize splif" %(splif_end-splif_begin))
      return({(0,0): splif_array})

    if "combined" in self.feature_types:
      begin = time.time()
      combined_array = self.concatenate_features(protein_xyz, protein_ob, ligand_xyz, ligand_ob)
      end = time.time()
      print("Took %f seconds to perform combined featurization" %(end-begin))
      return({(0,0): combined_array})

    if ligand_only:
      system = ligand 
      system_xyz = ligand_xyz
    else:
      system_xyz, system_ob = self.merge_molecules(protein_xyz, protein, ligand_xyz, ligand)
    print(time.time()-a)
    original_system = deepcopy(system)
    original_system.save_pdb(system_pdb_file)

    self.box = self.generate_box(original_system)
    print(time.time()-a)
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

    for key, transformed_system in transformed_systems.iteritems(): 
      transformed_box = self.generate_box(transformed_system)
      transformed_boxes[key] = transformed_box

    print(time.time()-a)

    if(save_intermediates):
      for key, transformed_system in transformed_systems.iteritems():
        print(key)
        print("%s/%s_%d_%d_system.pdb" %(save_dir, protein_name, key[0], key[1]))
        transformed_system.save_pdb("%s/%s_%d_%d_system.pdb" %(save_dir, protein_name, key[0], key[1]))
        pickle.dump(transformed_system, open("%s/%s_%d_%d_system.pickle" %(save_dir, protein_name, key[0], key[1]), "wb"))
      for key, transformed_box in transformed_boxes.iteritems():
        transformed_box.save_pdb("%s/%s_%d_%d_box.pdb" %(save_dir, protein_name, key[0], key[1]))
        pickle.dump(transformed_box, open("%s/%s_%d_%d_box.pickle" %(save_dir, protein_name, key[0], key[1]), "wb"))

    return(features)

  def get_xyz_from_ob(self, ob_mol):
    '''
    returns an m x 3 np array of 3d coords
    of given openbabel molecule
    '''

    xyz = np.zeros((ob_mol.NumAtoms(), 3))
    for i, atom in enumerate(ob.OBMolAtomIter(ob_mol)):
      xyz[i,0] = atom.x()
      xyz[i,1] = atom.y()
      xyz[i,2] = atom.z()
    return(xyz)

  def load_molecule(self, molecule_file, remove_hydrogens=True):
    '''
    given molecule_file, returns a tuple of xyz coords of molecule 
    and an openbabel object representing that molecule 
    '''

    if ".mol2" in molecule_file:
      obConversion = ob.OBConversion()
      obConversion.SetInAndOutFormats("mol2", "pdb")
      ob_mol = ob.OBMol()
      obConversion.ReadFile(ob_mol, molecule_file)
    else:
      obConversion = ob.OBConversion()
      obConversion.SetInAndOutFormats("pdb", "pdb")
      ob_mol = ob.OBMol()
      obConversion.ReadFile(ob_mol, molecule_file)

    #ob_mol.DeleteHydrogens()
    ob_mol.AddHydrogens()
    xyz = self.get_xyz_from_ob(ob_mol)

    return xyz, ob_mol

  def generate_box(self, mol):
    '''
    generate_box takes as input a molecule of class PDB and removes all atoms outside of the given box dims
    '''

    molecule = deepcopy(mol)
    atoms_to_keep = []
    all_atoms = [a for a in molecule.topology.atoms]
    for atom in all_atoms:
      coords = np.abs(molecule.xyz[0][atom.index,:])
      if coords[0] <= (self.box_x/2.) and coords[1] <= (self.box_y/2.) and coords[2] <= (self.box_z/2.):
        atoms_to_keep.append(atom.index)
    return(molecule.atom_slice(atoms_to_keep))

  def subtract_centroid(self, xyz, centroid):
    '''
    subtracts the centroid, a numpy array of dim 3, from all coordinates of all atoms in the molecule
    '''

    xyz -= np.transpose(centroid)
    return(xyz)

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

  def merge_molecules(self, protein_xyz, protein, ligand_xyz, ligand):
    '''
    Takes as input protein and ligand objects of class PDB and adds ligand atoms to the protein,
    and returns the new instance of class PDB called system that contains both sets of atoms.
    '''

    system_xyz = np.array(np.vstack(np.vstack((protein_xyz, ligand_xyz))))

    system_ob = ob.OBMol(protein_ob)
    system_ob += ligand_ob

    return system_xyz, system_ob

  '''
  Adapted from: http://baoilleach.blogspot.com/2008/02/calculate-circular-fingerprints-with.html
  '''
  def bfs(self, mol, startatom, D):
    '''
    given openbabel molecule and a starting atom of type OBAtom, 
    finds all bonds out to degree D via a breath-first search 
    '''

    visited_atoms = set()
    atoms_to_add = []
    bonds_to_add = []
    queue = deque([(startatom, 0)])
    while queue:
      atom, depth = queue.popleft()
      index = atom.GetIndex()
      visited_atoms.add(index)
      atomtype = atom.GetType()
      if depth < D:
        for bond in ob.OBAtomBondIter(atom):
          if bond not in bonds_to_add: bonds_to_add.append(bond)
      if depth < D:
        for atom in ob.OBAtomAtomIter(atom):
          if atom.GetIndex() not in visited_atoms:
            queue.append((atom, depth+1))
    return(bonds_to_add)

  def construct_fragment_from_bonds(self, bonds):
    '''
    takes as input a list of bonds of type OBBond and constructs a new 
    openbabel molecule from those bonds and the atoms that constitute 
    the start and end of those bonds. 
    '''

    fragment = ob.OBMol()
    added_atoms = []
    
    for bond in bonds: 
      atom_i = bond.GetBeginAtom()
      atom_j = bond.GetEndAtom()

      if atom_i not in added_atoms:
        fragment.AddAtom(atom_i)
        added_atoms.append(atom_i)
      atom_i_index = added_atoms.index(atom_i)

      if atom_j not in added_atoms:
        fragment.AddAtom(atom_j)
        added_atoms.append(atom_j)
      atom_j_index = added_atoms.index(atom_j)

      fragment.AddBond(atom_i_index + 1, atom_j_index + 1, bond.GetBondOrder())

    for i, fragment_bond in enumerate(ob.OBMolBondIter(fragment)):
      mol_bond = bonds[i]
      if mol_bond.IsAromatic(): fragment_bond.SetAromatic()

    return(fragment)

  def compute_ecfp(self, system_ob, start_atom, max_degree=2):
    '''
    Given an openbabel molecule and a starting atom (OBAtom object),
    compute the ECFP[max_degree]-like representation for that atom. 
    Returns (for now) a SMILES string representing the resulting fragment.

    TODO(enf): Optimize this! Try InChi key and other approaches
    to improving this representation. 
    '''

    fragment = ob.OBMol()
    
    bonds_to_add = self.bfs(system_ob, start_atom, max_degree)
    fragment = self.construct_fragment_from_bonds(bonds_to_add)
    obConversion = ob.OBConversion()
    obConversion.SetOutFormat("can")
    smiles = obConversion.WriteString(fragment).split("\t")[0]
    return(smiles)

  def hash_ecfp(self, ecfp, power):
    '''
    Returns an int of size 2^power representing that 
    ECFP fragment. Input must be a string. 
    '''

    md5 = hashlib.md5()
    md5.update(ecfp)
    digest = md5.hexdigest()
    ecfp_hash = int(digest, 16) % (2 ** power)
    return(ecfp_hash)

  def hash_ecfp_pair(self, ecfp_pair, power):
    '''
    Returns an int of size 2^power representing that 
    ECFP pair. Input must be a tuple of strings. 
    '''

    ecfp = "%s,%s" %(ecfp_pair[0],ecfp_pair[1])
    md5 = hashlib.md5()
    md5.update(ecfp)
    digest = md5.hexdigest()
    ecfp_hash = int(digest, 16) % (2 ** power)
    return(ecfp_hash)


  def compute_all_ecfp(self, system_ob, indices = None, degree=2):
    '''
    For each atom:
      Obtain molecular fragment for all atoms emanating outward to given degree. 
      For each fragment, compute SMILES string (for now) and hash to an int.
      Return a dictionary mapping atom index to hashed SMILES.
    '''

    ecfp_dict = {}

    for atom in ob.OBMolAtomIter(system_ob):
      if indices is not None:
        if atom.GetIndex() not in indices: continue
      ecfp_dict[atom.GetIndex()] = "%s,%s" %(atom.GetType(), self.compute_ecfp(system_ob, atom, degree))

    return(ecfp_dict)

  def compute_ecfp_features(self, system_ob, ecfp_degree, ecfp_power):
    '''
    Takes as input an openbabel molecule, ECFP radius, and number of bits to store 
    ECFP features (2^ecfp_power will be length of ECFP array);
    Returns an array of size 2^ecfp_power where array at index i has a 1 if that ECFP fragment
    is found in the molecule and array at index j has a 0 if ECFP fragment not in molecule. 
    '''

    ecfp_dict = self.compute_all_ecfp(system_ob, degree=ecfp_degree)
    ecfp_vec = [self.hash_ecfp(ecfp, self.ecfp_power) for index, ecfp in ecfp_dict.iteritems()]
    ecfp_array = np.zeros(2 ** self.ecfp_power)
    ecfp_array[sorted(ecfp_vec)] = 1.0
    return(ecfp_array)

  def featurize_binding_pocket_ecfp(self, protein_xyz, protein, ligand_xyz, ligand, pairwise_distances = None, cutoff = 4.5):
    '''
    Computes array of ECFP for both the ligand and the binding pocket region of the protein.
    '''

    if pairwise_distances is None: pairwise_distances = compute_pairwise_distances(protein_xyz, ligand_xyz)
    contacts = np.nonzero((pairwise_distances < 4.5))
    protein_atoms = set([int(c) for c in contacts[0].tolist()])
    protein_ecfp_dict = self.compute_all_ecfp(protein, indices = protein_atoms, degree = self.ecfp_degree)
    ligand_ecfp_dict = self.compute_all_ecfp(ligand, degree = self.ecfp_degree)

    protein_vec = [self.hash_ecfp(ecfp, self.ecfp_power) for index, ecfp in protein_ecfp_dict.iteritems()]
    ligand_vec = [self.hash_ecfp(ecfp, self.ecfp_power) for index, ecfp in ligand_ecfp_dict.iteritems()]

    protein_array = np.zeros(2 ** self.ecfp_power)
    ligand_array = np.zeros(2 ** self.ecfp_power)

    protein_array[sorted(protein_vec)] = 1.0
    print("Protein binding pocket has %d non-zero bits" %(len([f for f in protein_array if f == 1.])))
    ligand_array[sorted(ligand_vec)] = 1.0
    print("Ligand binding pocket has %d non-zero bits" %(len([f for f in ligand_array if f == 1.])))

    return(np.concatenate([ligand_array, protein_array]))



  def featurize_splif(self, protein_xyz, protein, ligand_xyz, ligand, pairwise_distances = None):
    '''
    Computes a SPLIF array of size 2^(self.splif_power) by first computing all protein atoms that come into contact
    with the ligand, then computing all ECFP fragments for all ligand atoms and all contacting protein atoms, 
    and then hashing to an int all contacting ECFP fragments. 
    '''

    if pairwise_distances is None: pairwise_distances = compute_pairwise_distances(protein_xyz, ligand_xyz)
    contact_bins = [(0,2.0), (2.0,3.0), (3.0,4.5)]
    splif_arrays = []
    protein_contacts = set()
    for contact_bin in contact_bins:
      contacts = np.nonzero((pairwise_distances > contact_bin[0]) & (pairwise_distances < contact_bin[1]))
      protein_atoms = set([int(c) for c in contacts[0].tolist()])
      protein_contacts = protein_contacts.union(protein_atoms)
      protein_ecfp_dict = self.compute_all_ecfp(protein, indices = protein_atoms, degree = self.ecfp_degree)
      ligand_ecfp_dict = self.compute_all_ecfp(ligand, degree = self.ecfp_degree)
      contacts = zip(contacts[0], contacts[1])
      ecfp_pairs = [(protein_ecfp_dict[contact[0]], ligand_ecfp_dict[contact[1]]) for contact in contacts]
      splif_vec = [self.hash_ecfp_pair(ecfp_pair, self.splif_power) for ecfp_pair in ecfp_pairs]
      splif_array = np.zeros(2 ** self.splif_power)
      splif_array[sorted(splif_vec)] = 1
      splif_arrays.append(splif_array)
      print("Splif array bin has %d non-zero entries" %(len([f for f in splif_array if f == 1.])))
    print("there are %d protein_atoms to ecfp" %(len(protein_contacts)))

    splif = np.concatenate(splif_arrays)
    return(splif)

  def is_hydrogen_bond(self, protein_xyz, protein, ligand_xyz, ligand, contact):
    '''
    Determine if a pair of atoms (contact = tuple of protein_atom_index, ligand_atom_index)
    between protein and ligand represents a hydrogen bond. Returns a boolean result. 
    '''

    protein_atom_index = contact[0]
    ligand_atom_index = contact[1]
    protein_atom = protein.GetAtomById(protein_atom_index)
    ligand_atom = ligand.GetAtomById(ligand_atom_index)
    if protein_atom.IsHbondAcceptor() and ligand_atom.IsHbondDonor():
      for atom in ob.OBAtomAtomIter(ligand_atom):
        if atom.GetAtomicNum() == 1:
          hydrogen_xyz = ligand_xyz[atom.GetIndex(),:]
          v1 = protein_xyz[protein_atom_index,:] - hydrogen_xyz
          v2 = ligand_xyz[ligand_atom_index,:] - hydrogen_xyz
          angle = angle_between(v1,v2) * 180. / np.pi
          if(angle > (180 - self.hbond_angle_cutoff) and angle < (180. + self.hbond_angle_cutoff)):
            return True
    elif ligand_atom.IsHbondAcceptor() and protein_atom.IsHbondDonor():
      for atom in ob.OBAtomAtomIter(protein_atom):
        if atom.GetAtomicNum() == 1:
          hydrogen_xyz = protein_xyz[atom.GetIndex(),:]
          v1 = protein_xyz[protein_atom_index,:] - hydrogen_xyz
          v2 = ligand_xyz[ligand_atom_index,:] - hydrogen_xyz
          angle = angle_between(v1,v2) * 180. / np.pi
          if(angle > (180 - self.hbond_angle_cutoff) and angle < (180. + self.hbond_angle_cutoff)):
            return True
    return False     

  def featurize_hydrogen_bond_array(self, protein_xyz, protein, ligand_xyz, ligand, pairwise_distances = None):
    '''
    Given protein and ligand, determine all hydrogen bonds between them. For each ligand_ecfp_i, protein_ecfp_j pair,
    determine if there is a hydrogen bond between them, and, if so, set to 1 the entry in an array corresponding 
    to the integer index of the hash of that ecfp pair. Return resulting array of size 2^self.splif_power.

    TODO(enf): Change this to a 3-bit vector, where each bit is number of hydrogen bonds for each bin. These
    bins will be defined as in: 
    http://evans.rc.fas.harvard.edu/pdf/smnr_2009_Kwan_Eugene.pdf
    '''

    if pairwise_distances is None: pairwise_distances = compute_pairwise_distances(protein_xyz, ligand_xyz)

    contacts = np.nonzero(pairwise_distances < self.hbond_dist_cutoff)
    protein_atoms = set([int(c) for c in contacts[0].tolist()])
    protein_ecfp_dict = self.compute_all_ecfp(protein, indices = protein_atoms, degree = self.ecfp_degree)
    ligand_ecfp_dict = self.compute_all_ecfp(ligand, degree = self.ecfp_degree)
    contacts = zip(contacts[0], contacts[1])
    hydrogen_bond_contacts = []
    for contact in contacts:
      if self.is_hydrogen_bond(protein_xyz, protein, ligand_xyz, ligand, contact):
        hydrogen_bond_contacts.append(contact)

    ecfp_pairs = [(protein_ecfp_dict[contact[0]], ligand_ecfp_dict[contact[1]]) for contact in hydrogen_bond_contacts]
    splif_vec = [self.hash_ecfp_pair(ecfp_pair, self.splif_power) for ecfp_pair in ecfp_pairs]
    splif_array = np.zeros(2 ** self.splif_power)
    splif_array[sorted(splif_vec)] = 1.
    print("There are %d h-bonds" %(len([f for f in splif_array if f == 1.])))
    return(splif_array)

  #TODO(enf): add flag to specify which features to perform and then concat.
  def concatenate_features(self, protein_xyz, protein, ligand_xyz, ligand):
    '''
    Compute all feature types and concatenate them into one vector. 
    '''

    a = time.time()
    pairwise_distances = compute_pairwise_distances(protein_xyz, ligand_xyz)
    print("Computing pairwise_distances took %f seconds" %(time.time()-a))
    a = time.time()
    binding_pocket_ecfp_features = self.featurize_binding_pocket_ecfp(protein_xyz, protein, ligand_xyz, ligand, pairwise_distances = pairwise_distances)
    print("Computing binding_pocket_ecfp_features took %f seconds" %(time.time()-a))
    a = time.time()
    splif_features = self.featurize_splif(protein_xyz, protein, ligand_xyz, ligand, pairwise_distances)
    print("computing splif_features took %f seconds" %(time.time()-a))
    a = time.time()
    #hydrogen_bond_features = self.featurize_hydrogen_bond_array(protein_xyz, protein, ligand_xyz, ligand, pairwise_distances = pairwise_distances)
    #print("Computing hydrogen bonds took %f seconds" %(time.time()-a))

    features = np.concatenate([binding_pocket_ecfp_features, splif_features])
    return(features)