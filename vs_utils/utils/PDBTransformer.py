from nnscore_pdb import PDB 
from copy import deepcopy
import numpy as np
import pickle 

'''Pseudocode 

1. Create objects of class PDB for protein and ligand
2. Compute centroid of ligand
3. translate coordinates of all other atoms such that ligand centroid is now origin, for both ligand PDB and protein PDB
4. combine PDB objects
5. Write out new PDB 
'''

'''
Pseudocode, detail: 

1. Create objects of class PDB for protein and ligand 
2. Compute centroid of ligand
	a. define a function for computing centroid given a PDB object. Maybe extend to make it accept just coordinates?
		-further decomposition:
			Take object as input --> convert to coordinates --> then send to separate method that computes centroid 
3. For each atom for each PDB oject: 
	-change x, y, and z coordinates by subtracting out centroid x, y, and z
4. Function that combines PDB objects: 
	def mergePDBs(pdb1, pdb2)
		copy pdb1
		add all from pdb2 to pdb1, but offset atom indices by max(atom index of pdb 1)
'''

def compute_centroid(molecule):
	coordinates = [molecule.all_atoms[atom].coordinates.coords for atom in molecule.all_atoms]
	centroid = np.sum(coordinates, axis=0)/float(len(coordinates))
	return(centroid)

class PDBTransformer:
	def __init__(self):
		self.protein = None
		self.ligand = None
		self.system = None
		self.box = None
		self.box_x = 100000.0
		self.box_y = 100000.0
		self.box_z = 100000.0


	def transform(self, protein_pdb, protein_pdbqt, ligand_pdb, ligand_pdbqt, system_pdb, box_pdb, box_pickle, box_x, box_y, box_z):
		self.protein = PDB()
		self.protein.load_from_files(protein_pdb, protein_pdbqt)

		self.ligand = PDB()
		self.ligand.load_from_files(ligand_pdb, ligand_pdbqt)

		self.box_x = box_x
		self.box_y = box_y
		self.box_z = box_z

		ligand_centroid = compute_centroid(self.ligand)
		self.ligand = self.subtract_centroid(self.ligand, ligand_centroid)
		print(self.ligand)
		print(self.ligand.all_atoms[1])
		self.protein = self.subtract_centroid(self.protein, ligand_centroid)

		self.system = self.merge_molecules(self.protein, self.ligand, self.system)
		print(self.system)
		print(self.system.all_atoms[1])
		#for atomindex in self.system.all_atoms:
		#	print(self.system.all_atoms[atomindex])
		#	print(self.system.all_atoms[atomindex].create_pdb_line(atomindex))
		self.system.save_pdb(system_pdb)

		self.box = self.generate_box(self.system)
		self.box.save_pdb(box_pdb)
		pickle.dump(self.box, open(box_pickle, "wb"))


	def generate_box(self, molecule):
		atoms_to_remove = []
		for index, atom in molecule.all_atoms.iteritems():
			coords = np.abs(atom.coordinates.coords)
			if coords[0] >= (self.box_x/2.) or coords[1] >= (self.box_y/2.) or coords[2] >= (self.box_z/2.):
				atoms_to_remove.append((index, atom))

		for index, atom in atoms_to_remove:
			molecule = self.remove_atom(molecule, index, atom)

		return(molecule)

	def remove_atom(self, molecule, atom_index, atom):
		if molecule.all_atoms[atom_index] != atom:
			print("Removing atoms from dictionary not safe!!")
			return

		del molecule.all_atoms[atom_index]
		return(molecule)

	def subtract_centroid(self, molecule, centroid):
		for atom in molecule.all_atoms:
			molecule.all_atoms[atom].coordinates.coords -= centroid
		return(molecule)

	def merge_molecules(self, protein, ligand, system):
		system = deepcopy(protein)
		greatest_index = len(protein.all_atoms) + 1
		autoindex = greatest_index + 1
		#print(ligand.all_atoms)
		for index, ligand_atom in ligand.all_atoms.iteritems():
			#print(ligand_atom)
			#system.all_atoms[autoindex] = ligand_atom
			#autoindex += 1
			system.add_new_atom(ligand_atom)
		#system.add_new_atoms(ligand_atoms)
		print system.all_atoms[greatest_index+1]
		return(system)


