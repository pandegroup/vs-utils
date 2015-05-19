"""
Prep PDBBind molecules for processing by nnscore.
"""
import argparse
import openbabel
import os
from vs_utils.scripts.public_data.nnscore_utilities import hydrogenate_and_compute_partial_charges

def parse_args(input_args=None):
  """Parse command-line arguments."""
  parser = argparse.ArgumentParser()
  parser.add_argument('--pdbbind-dir', required=1,
                      help='Directory containing pdbbind data.')
  return parser.parse_args(input_args)

def main(pdbbind_dir):
  assert os.path.isdir(pdbbind_dir)
  subdirs = [d for d in os.listdir(pdbbind_dir) if
      os.path.isdir(os.path.join(pdbbind_dir, d))]
  print subdirs
  for d in subdirs:
    print "processing %s" % d
    subdir = os.path.join(pdbbind_dir, d)
    ligand, protein = None, None
    for f in os.listdir(subdir):
      if "_ligand.mol2" in f:
        print "ligand_name: %s" % f
        ligand = f
      elif "_protein.pdb" in f:
        print "protein_name: %s" %f
        protein = f
    if not ligand or not protein:
      raise ValueError("Ligand or Protein missing in %s" % d)
    ligand_file = os.path.join(subdir, ligand)
    protein_file = os.path.join(subdir, protein)
    hydrogenate_and_compute_partial_charges(ligand_file, "mol2", subdir)
    hydrogenate_and_compute_partial_charges(protein_file, "pdb", subdir)

if __name__ == '__main__':
  args = parse_args()
  main(args.pdbbind_dir)
