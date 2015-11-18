"""
Featurize PDBBind.
"""
import argparse
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

def featurize_pdbbind(pdbbind_dir, pickle_out):
  """Featurize all entries in pdbbind_dir and write features to pickle_out

  pdbbind_dir should be a dir, with K subdirs, one for each protein-ligand
  complex to be featurized. The ligand and receptor should each have a pdb
  and pdbqt file. The ligand files should end in '_ligand_hyd.${FILETYPE}'
  while the receptor files should end in '_protein_hyd.${FILETYPE}'

  pdbbind_dir: string
    Path to pdbbind directory.
  pickle_out: string
    Path to write pickle output.
  """
  assert os.path.isdir(pdbbind_dir)
  # Instantiate copy of binana vector
  binana = Binana()
  feature_vectors = {}

  # Extract the subdirectories in pdbbind_dir
  subdirs = [d for d in os.listdir(pdbbind_dir) if
      os.path.isdir(os.path.join(pdbbind_dir, d))]
  # TODO(rbharath): ONLY FOR DEBUGGING!

  num_atoms = len(Binana.atom_types)
  # See features/tests/nnscore_test.py:TestBinana.testComputeInputVector
  # for derivation.
  feature_len = (3*num_atoms*(num_atoms+1)/2 + num_atoms + 12 + 6 + 3 + 6 +
      3 + 6 + 3 + 1)
  for count, d in enumerate(subdirs):
    print "\nprocessing %d-th pdb %s" % (count, d)
    subdir = os.path.join(pdbbind_dir, d)

    print "About to extract ligand and protein input files"
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

    print "Extracted Input Files:"
    print (ligand_pdb, ligand_pdbqt, protein_pdb, protein_pdbqt)
    if (not ligand_pdb or not ligand_pdbqt or not protein_pdb or not
        protein_pdbqt):
        raise ValueError("Required files not present for %s" % d)
    ligand_pdb_path = os.path.join(subdir, ligand_pdb)
    ligand_pdbqt_path = os.path.join(subdir, ligand_pdbqt)
    protein_pdb_path = os.path.join(subdir, protein_pdb)
    protein_pdbqt_path = os.path.join(subdir, protein_pdbqt)

    print "About to load ligand from input files"
    ligand_pdb_obj = PDB()
    ligand_pdb_obj.load_from_files(ligand_pdb_path, ligand_pdbqt_path)

    print "About to load protein from input files"
    protein_pdb_obj = PDB()
    protein_pdb_obj.load_from_files(protein_pdb_path, protein_pdbqt_path)

    print "About to generate feature vector."
    vector = binana.compute_input_vector(ligand_pdb_obj,
        protein_pdb_obj)
    feature_vectors[d] = vector
    if len(vector) != feature_len:
      raise ValueError("Feature length incorrect on %s" % d)
    print "Feature vector generated correctly."

  with open(pickle_out, "wb") as f:
    pickle.dump(feature_vectors, f)

if __name__ == '__main__':
  args = parse_args()
  featurize_pdbbind(args.pdbbind_dir, args.pickle_out)
