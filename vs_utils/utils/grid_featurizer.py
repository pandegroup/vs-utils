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

  coordinates = molecule.xyz 
  centroid = np.mean(coordinates, axis=1)
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


class grid_featurizer:
  def __init__(self):
    self.protein = None
    self.protein_graph = None 
    self.protein_ob = None 
    self.ligand = None
    self.ligand_graph = None
    self.ligand_ob = None
    self.system = None
    self.system_graph = None
    self.system_ob = None
    self.box = None
    self.ecfp_power = 10
    self.ecfp_degree = 2

    self.box_x = 100000.0
    self.box_y = 100000.0
    self.box_z = 100000.0


  def transform(self, protein_pdb, ligand_pdb, save_dir, box_x=16.0, box_y=16.0, box_z=16.0, 
    nb_rotations = 0, nb_reflections=0, feature_types="ecfp", 
    ecfp_degree=2, ecfp_power = 10, save_intermediates=False, ligand_only = False, **kwargs):
    '''Takes as input files (strings) for pdb/pdbqt of the protein, pdb/pdbqt of the ligand, a filename to be saved 
    of the merged system called system_pdb, a filename to be saved of the box called box_pdb, a filename with a .p extension
    of the box in pickle format called box_pickle, and 3 floats for the x,y,z dimensions of the box in Angstroms

    This function then computes the centroid of the ligand; decrements this centroid from the atomic coordinates of protein and
    ligand atoms, and then merges the translated protein and ligand. This combined system/complex is then saved.

    This function then removes all atoms outside of the given box dimensions, and saves the resulting box in PDB file format
    as well as in a pickle format of the box's instance of the PDB object.
    '''
    a = time.time()
    protein_name = str(protein_pdb).split("/")[len(str(protein_pdb).split("/"))-2]
    system_pdb_file = "%s/%s.pdb" %(save_dir, protein_name)
    system_pickle_file = "%s/%s.pickle" %(save_dir, protein_name)

    features = []

    self.box_x = float(box_x)/10.0
    self.box_y = float(box_y)/10.0
    self.box_z = float(box_z)/10.0

    self.ecfp_degree = ecfp_degree
    self.ecfp_power = ecfp_power
    print("Starting")
    if not ligand_only:
      self.protein, self.protein_ob = self.load_molecule(protein_pdb) 
    print("Loading ligand")
    self.ligand, self.ligand_ob = self.load_molecule(ligand_pdb)
    print("Loaded ligand")

    if "ecfp" in feature_types:
      ecfp_array = self.compute_ecfp_features(self.ligand_ob, ecfp_degree, ecfp_power)
      with gzip.open("%s/%s_ecfp_array.pkl.gz" %(save_dir, protein_name), "wb") as f:
        pickle.dump(ecfp_array, f)
      return({(0,0): ecfp_array})

    centroid = compute_centroid(self.ligand)
    self.ligand = self.subtract_centroid(self.ligand, centroid)
    if not ligand_only:
      self.protein = self.subtract_centroid(self.protein, centroid)
    print(time.time()-a)
    if ligand_only:
      self.system = self.ligand 
      self.system_ob = self.ligand_ob
    else:
      self.system, self.system_ob = self.merge_molecules(self.protein, self.protein_ob, self.ligand, self.ligand_ob)

    print(time.time()-a)
    original_system = deepcopy(self.system)
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



      #for key, transformed_system in transformed_systems.iteritems():
      #  ecfp_dict[key] = self.compute_ecfp(transformed_system)

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

  def load_molecule(self, molecule_file):
    if ".pdb" in molecule_file:
      molecule = md.load(molecule_file)
    else:
      molecule = None
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
    
    '''
    G = nx.Graph()

    for bond in ob.OBMolBondIterator(ob_mol):
      atom0 = bond.GetBeginAtom()
      atom1 = bond.GetEndAtom()

      G.add_edge(atom0.GetIndex(), atom1.GetIndex())
      if bond.IsAromatic():
        G[atom0.GetIndex()][atom1.GetIndex()]["BondOrder"] = 1.5
      else:
        G[atom0.GetIndex()][atom1.GetIndex()]["BondOrder"] = bond.GetBondOrder()

      G.node[atom0.GetIndex()]["AtomicNumber"] = atom0.GetAtomicNum()
      G.node[atom1.GetIndex()]["AtomicNumber"] = atom1.GetAtomicNum()
      G.node[atom0.GetIndex()]["SybylType"] = atom0.GetType()
      G.node[atom1.GetIndex()]["SybylType"] = atom1.GetType()
    '''

    return molecule, ob_mol




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
    molecule.xyz -= centroid 
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

  def merge_molecules(self, protein, protein_ob, ligand, ligand_ob):
    '''
    Takes as input protein and ligand objects of class PDB and adds ligand atoms to the protein,
    and returns the new instance of class PDB called system that contains both sets of atoms.
    '''

    system = deepcopy(protein)
    ligand_atoms = [a for a in ligand.topology.atoms]
    ligand_residue = [r for r in ligand.topology.residues][0]
    ligand_residue_name = ligand_residue.name 
    ligand_residue_resSeq = ligand_residue.resSeq 
    new_chain = system.topology.add_chain()
    new_residue = system.topology.add_residue(ligand_residue_name, chain = new_chain, resSeq = ligand_residue_resSeq) 
    max_serial = ligand_atoms[len(ligand_atoms)-1].serial
    for i, atom in enumerate(ligand_atoms):
      system.topology.add_atom(atom.name, atom.element, residue = new_residue, serial = max_serial + i + 1)

    system.xyz = np.array(np.vstack(np.vstack((system.xyz[0], ligand.xyz[0]))))
    '''
    num_protein_atoms = len([a for a in protein.topology.atoms])
    mapping = {i: (i+num_protein_atoms) for i in range(0, len(ligand_atoms))}
    ligand_graph=nx.relabel_nodes(ligand_graph,mapping)

    system_graph = nx.union(protein_graph, ligand_graph)
    '''
    system_ob = ob.OBMol(protein_ob)
    system_ob += ligand_ob

    return system, system_ob

  '''
  Adapted from: http://baoilleach.blogspot.com/2008/02/calculate-circular-fingerprints-with.html
  '''
  def bfs(self, mol, startatom, D):
    visited = [False] * (mol.NumAtoms())
    atoms_to_add = []
    bonds_to_add = []
    queue = deque([(startatom, 0)])
    while queue:
        atom, depth = queue.popleft()
        index = atom.GetIndex()
        visited[index] = True
        atomtype = atom.GetType()
        if depth < D:
          for bond in ob.OBAtomBondIter(atom):
            if bond not in bonds_to_add: bonds_to_add.append(bond)
        if depth < D:
            for atom in ob.OBAtomAtomIter(atom):
               if not visited[atom.GetIndex()]:
                   queue.append((atom, depth+1))
    return(bonds_to_add)

  def construct_fragment_from_bonds(self, bonds):
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
    fragment = ob.OBMol()
    
    bonds_to_add = self.bfs(system_ob, start_atom, max_degree)
    fragment = self.construct_fragment_from_bonds(bonds_to_add)
    obConversion = ob.OBConversion()
    obConversion.SetOutFormat("smiles")
    smiles = obConversion.WriteString(fragment)
    return(smiles)

  def hash_ecfp(self, ecfp, power):
    #ecfp_hash = abs(hash(ecfp)) % (2 ** power)
    md5 = hashlib.md5()
    md5.update(ecfp)
    digest = md5.hexdigest()
    ecfp_hash = int(digest, 16) % (2 ** power)
    return(ecfp_hash)

  def compute_all_ecfp(self, system_ob, degree=2):
    '''
    For each atom:
      Obtain bondgraph for all atoms emanating outward to specified degree.
      Save that bondgraph in a dict w/ key as atom index (is this safe?) and val as bondgr
    '''
    ecfp_dict = {}
    for atom in ob.OBMolAtomIter(system_ob):
      ecfp_dict[atom.GetIndex()] = self.compute_ecfp(system_ob, atom, degree)

    return(ecfp_dict)


  def compute_ecfp_features(self, system_ob, ecfp_degree, ecfp_power):
    ecfp_dict = self.compute_all_ecfp(system_ob, ecfp_degree)
    #print([fp for idx, fp in ecfp_dict.iteritems()][0:20])
    ecfp_vec = [self.hash_ecfp(ecfp, self.ecfp_power) for index, ecfp in ecfp_dict.iteritems()]
    ecfp_array = np.zeros(2 ** self.ecfp_power)
    ecfp_array[sorted(ecfp_vec)] = 1.0
    return(ecfp_array)

  def featurize_ecfp(self, system, degree=2):
    '''


    '''


