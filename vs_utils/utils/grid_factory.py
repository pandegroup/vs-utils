from nnscore_pdb import PDB 
from copy import deepcopy
import numpy as np
import pickle
import csv
import os
from collections import defaultdict

class GridGenerator:
  def __init__(self):
    self.box_x = 16.0
    self.box_y = 16.0
    self.box_z = 16.0
    self.x_linspace = None
    self.y_linspace = None
    self.z_linspace = None
    self.box = None
    self.grid_res = 0.5
    self.grid = {}
    self.grid_vector = []
    self.num_features = 3

    path = os.path.dirname(os.path.realpath(__file__))

    with open(os.path.join(path, 'info/periodic_table.p'), 'rb') as handle:
      self.periodic_table = pickle.load(handle)
    self.periodic_table = dict((k.upper(),v) for k,v in self.periodic_table.iteritems())

    with open(os.path.join(path, 'info/protein_heavy_atoms.p'), 'rb') as handle:
      self.protein_heavy_atoms = pickle.load(handle)

    with open(os.path.join(path, 'info/all_atom_names.p'), 'rb') as handle:
      self.all_atom_names = pickle.load(handle)
    self.terminal_atoms = ["H1", "1H", "HN1", "HN1","HT1","HNCAP","H2", "2H", "HN2","HN2","HT2","H3","3H","HN3","HN3","HT3",
                "OT1", "OXT", "OT2"]
    self.all_atom_names += self.terminal_atoms
    self.all_atom_names = list(set(self.all_atom_names))


  def transform(self, box, box_x, box_y, box_z, voxel_width, grid_pickle_file, num_features=3):
    '''
    Takes as input a pickle file of a "box" generated from PDBTransformer's transform() function. This "box" will be an
    instance of class PDB defined in nnscore_pdb.py. It will contain all atoms within a rectangular prism with edge lengths 
    defined as box_x, box_y, and box_z. These dimensions must be the same as the dimensions fed to PDBTransformer.
    [To do: maybe make "box" a class that inherits PDB?]

    You must also provide: voxel_width, a float that describes in angstroms the edge length of a single voxel.
    grid_pickle_file: pickle file for saving 4d np array of featurized grid
    num_features: this will specify the length of the 4th dimension of the tensor. Default is 3.

    '''

    self.box_x = float(box_x)
    self.box_y = float(box_y)
    self.box_z = float(box_z)
    self.grid_res = float(voxel_width)
    self.num_features = 3
    self.box = box


    self.initialize_grid()

    self.featurize_atoms()
    grid_valid = self.check_grid()
    if grid_valid:
      pickle.dump(self.grid, open(grid_pickle_file, "wb"))
    else:
      print("Something went wrong...")

    return(self.grid)



  def read_box_from_pdb(self, pdb_file):
    '''
    Functionality isn't supported yet. For now, you must pass in to transform() a pickle file that contains 
    an instance of class PDB described in nnscore_pdb that contains the box that will be converted to a grid. 
    '''
    return

  def initialize_grid(self):
    '''
    Initializes linspaces describing delimiters of voxels, and then initializes 4-tensor of featurized grid with zeros.
    '''
    self.x_linspace = np.linspace(-1.0 * self.box_x / 2.0, self.box_x / 2.0, (self.box_x + self.grid_res) / self.grid_res)
    self.y_linspace = np.linspace(-1.0 * self.box_y / 2.0, self.box_y / 2.0, (self.box_y + self.grid_res) / self.grid_res)
    self.z_linspace = np.linspace(-1.0 * self.box_z / 2.0, self.box_z / 2.0, (self.box_z + self.grid_res) / self.grid_res)

    self.grid = np.zeros((len(self.x_linspace)-1, len(self.y_linspace)-1, len(self.z_linspace)-1, self.num_features))

    return

  def featurize_atoms(self):
    for index, atom in self.box.all_atoms.iteritems():
      atomname = atom.atomname.strip()
      atomtype = atom.atomtype.strip()
      resname = atom.residue.strip()

      if resname in self.box.protein_resnames:
        res_int = self.box.protein_resnames.index(resname) + 1
        try:
          atom_type_int = self.protein_heavy_atoms.index(atomname) + 1
          #atom_type_int = self.all_atom_names.index(atomname) + 1
        except:
          atom_type_int = 0
      else:
        res_int = 0
        atom_type_int = 0

      element = atom.element
      element_int = self.periodic_table[element]

      grid_index = self.locate_atom(atom.coordinates.coords)
      atom_features = [element_int, res_int, atom_type_int]
      if np.linalg.norm(self.grid[grid_index[0]][grid_index[1]][grid_index[2]], ord=0) > 0: print("WARNING: Atom already in this voxel!!")

      self.grid[grid_index[0]][grid_index[1]][grid_index[2]] = np.array(atom_features)
    print("Finished featurizing atoms")
    return

  def locate_atom(self, coords):
    
    for i in range(0, len(self.x_linspace)-1):
      if coords[0] >= self.x_linspace[i] and coords[0] <= self.x_linspace[i+1]: break
    for j in range(0, len(self.y_linspace)-1):
      if coords[1] >= self.y_linspace[j] and coords[1] <= self.y_linspace[j+1]: break
    for k in range(0, len(self.z_linspace)-1):
      if coords[2] >= self.z_linspace[k] and coords[2] <= self.z_linspace[k+1]: break

    location = [i, j, k]
    return(location)

  def check_grid(self):
    grid_non_none_features = 0
    target_non_none_features = len(self.box.all_atoms)
    for i in range(0, np.shape(self.grid)[0]):
      for j in range(0, np.shape(self.grid)[1]):
        for k in range(0, np.shape(self.grid)[2]):
          if np.linalg.norm(self.grid[i][j][k], ord = 0) > 0: grid_non_none_features += 1
          
    if grid_non_none_features == target_non_none_features: 
      print("Correct number of features in grid!")
      return True
    else:
      return False



