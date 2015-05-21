"""
Featurize the core set of PDBBind.
"""
import argparse
import openbabel
import os
import re
import cPickle as pickle
from vs_utils.features.nnscore import Binana
from vs_utils.utils.nnscore_pdb import PDB

def parse_args(input_args=None):
  """Parse command-line arguments."""
  parser = argparse.ArgumentParser()
  parser.add_argument('--pdbbind-dir', required=1,
                      help='Directory containing pdbbind data.')
  parser.add_argument('--pickle-out', required=1,
                      help='Path to output pickled featured vectors.')
  return parser.parse_args(input_args)

def main(pdbbind_dir, pickle_out):
  # Instantiate copy of binana vector
  binana = Binana()
  feature_vectors = {}

  assert os.path.isdir(pdbbind_dir)
  subdirs = [d for d in os.listdir(pdbbind_dir) if
      os.path.isdir(os.path.join(pdbbind_dir, d))]
  #subdirs = ['1r5y']

  N = len(Binana.atom_types)
  # See features/tests/nnscore_test.py:TestBinana.testComputeInputVector
  # for derivation.
  feature_len = 3*N*(N+1)/2 + N + 12 + 6 + 3 + 6 + 3 + 6 + 3 + 1
  for d in subdirs:
    print "\nprocessing %s" % d
    subdir = os.path.join(pdbbind_dir, d)
    ligand_pdb, ligand_pdbqt = None, None
    protein_pdb, protein_pdbqt = None, None
    for f in os.listdir(subdir):
      if re.search("_ligand_hyd.pdb$", f):
        ligand_pdb = f
      elif re.search("_ligand_hyd.pdbqt$", f):
        ligand_pdbqt = f
      elif re.search("_protein_hyd.pdb$", f):
        protein_pdb = f
      elif re.search("_protein_hyd.pdbqt$", f):
        protein_pdbqt = f
    print (ligand_pdb, ligand_pdbqt, protein_pdb, protein_pdbqt)
    if (not ligand_pdb or not ligand_pdbqt or not protein_pdb or not
        protein_pdbqt):
        raise ValueError("Required files not present for %s" % d)
    ligand_pdb_path = os.path.join(subdir, ligand_pdb)
    ligand_pdbqt_path = os.path.join(subdir, ligand_pdbqt)
    protein_pdb_path = os.path.join(subdir, protein_pdb)
    protein_pdbqt_path = os.path.join(subdir, protein_pdbqt)

    print "About to load ligand from pdbs"
    ligand_pdb_obj = PDB()
    ligand_pdb_obj.load_from_files(ligand_pdb_path, ligand_pdbqt_path)

    print "About to load protein from pdbs"
    protein_pdb_obj = PDB()
    protein_pdb_obj.load_from_files(protein_pdb_path, protein_pdbqt_path)

    print "About to generate feature vector"
    vector = binana.compute_input_vector(ligand_pdb_obj,
        protein_pdb_obj)
    feature_vectors[d] = vector
    if len(vector) != feature_len:
      raise ValueError("Feature length incorrect on %s" % d)

  with open(pickle_out, "wb") as f:
    pickle.dump(feature_vectors, f)

if __name__ == '__main__':
  args = parse_args()
  main(args.pdbbind_dir, args.pickle_out)
