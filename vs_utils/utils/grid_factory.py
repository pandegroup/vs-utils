from nnscore_pdb import PDB 
from copy import deepcopy
import numpy as np
import pickle

class GridGenerator:
	def __init__(self):
		self.box_x = 20.0
		self.box_y = 20.0
		self.box_z = 20.0
		self.x_linspace = None
		self.y_linspace = None
		self.z_linspace = None
		self.box = None
		self.grid_res = 1.0
		self.grid = {}
		self.grid_vector = []
		self.num_features = 3


	def read_box_from_pdb(self, pdb_file):
		return

	def transform(self, box_pickle, box_x, box_y, box_z, grid_res, num_features, grid_3d_pickle, grid_1d_pickle):
		self.box_x = box_x
		self.box_y = box_y
		self.box_z = box_z
		self.grid_res = grid_res
		self.num_features = 3
		self.box = pickle.load(open(box_pickle, "rb"))


		self.initialize_grid()

		self.featurize_atoms()
		self.vectorize_grid()
		grid_valid = self.check_grid()
		if grid_valid:
			pickle.dump(self.grid, open(grid_3d_pickle, "wb"))
			pickle.dump(self.grid_vector, open(grid_1d_pickle, "wb"))
		else:
			print("Something went wrong...")



	def initialize_grid(self):
		self.x_linspace = np.linspace(-1.0 * self.box_x / 2.0, self.box_x / 2.0, (self.box_x + self.grid_res) / self.grid_res)
		self.y_linspace = np.linspace(-1.0 * self.box_y / 2.0, self.box_y / 2.0, (self.box_y + self.grid_res) / self.grid_res)
		self.z_linspace = np.linspace(-1.0 * self.box_z / 2.0, self.box_z / 2.0, (self.box_z + self.grid_res) / self.grid_res)


		for i in range(0, len(self.x_linspace)-1):
			self.grid[i] = {}
			for j in range(0, len(self.y_linspace)-1):
				self.grid[i][j] = {}
				for k in range(0, len(self.z_linspace)-1):
					self.grid[i][j][k] = [None for f in range(0,self.num_features)]
		print("Finished initializing grid")
		return

	def featurize_atoms(self):
		for index, atom in self.box.all_atoms.iteritems():
			atomname = atom.atomname.strip()
			atomtype = atom.atomtype.strip()
			resname = atom.residue.strip()
			if resname not in self.box.protein_resnames:
				resname = "LIG"
				atomname = None
			grid_index = self.locate_atom(atom.coordinates.coords)
			print(grid_index)
			atom_features = [atomname, atomtype, resname]
			print(atom_features)
			if None not in self.grid[grid_index[0]][grid_index[1]][grid_index[2]]: print("Atom already in this voxel!!")
			self.grid[grid_index[0]][grid_index[1]][grid_index[2]] = atom_features
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

	def vectorize_grid(self):
		for i in range(len(self.grid)):
			for j in range(len(self.grid[i])):
				for k in range(len(self.grid[i][j])):
					self.grid_vector.append(self.grid[i][j][k])
		print("Finished vectorizing grid")
		return

	def check_grid(self):
		grid_non_none_features = 0
		target_non_none_features = len(self.box.all_atoms) * self.num_features
		for feature in self.grid_vector:
			if None not in feature: grid_non_none_features += len(feature)
		print(grid_non_none_features)
		print(target_non_none_features)
		if grid_non_none_features == target_non_none_features: 
			print("Correct number of features in grid!")
			return True
		else:
			return False



