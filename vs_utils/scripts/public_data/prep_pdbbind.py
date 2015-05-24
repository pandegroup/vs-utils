"""
Prep PDBBind molecules for processing by nnscore.
"""
import argparse
import os
from vs_utils.scripts.public_data.nnscore_utilities import hydrogenate_and_compute_partial_charges

def parse_args(input_args=None):
  """Parse command-line arguments."""
  parser = argparse.ArgumentParser()
  parser.add_argument('--pdbbind-dir', required=1,
                      help='Directory containing pdbbind data.')
  return parser.parse_args(input_args)

def preprocess_pdbbind(pdbbind_dir):
  """Preprocess pdbbind files for Binana."""
  assert os.path.isdir(pdbbind_dir)

  # Extract the subdirectories in pdbbind_dir
  subdirs = [d for d in os.listdir(pdbbind_dir) if
      os.path.isdir(os.path.join(pdbbind_dir, d))]

  print "About to preprocess following subdirectories:"
  print subdirs

  for count, d in enumerate(subdirs):
    print "Processing %d-th entry %s" % (count, d)
    subdir = os.path.join(pdbbind_dir, d)
    ligand, protein = None, None
    for f in os.listdir(subdir):
      if "_ligand.mol2" in f:
        print "Input ligand: %s" % f
        ligand = f
      elif "_protein.pdb" in f:
        print "Input protein: %s" % f
        protein = f
    if not ligand or not protein:
      raise ValueError("Ligand or Protein missing in %s" % d)
    ligand_file = os.path.join(subdir, ligand)
    protein_file = os.path.join(subdir, protein)

    print "About to preprocess ligand."
    hydrogenate_and_compute_partial_charges(ligand_file, "mol2", subdir)

    print "About to preprocess protein."
    hydrogenate_and_compute_partial_charges(protein_file, "pdb", subdir)

if __name__ == '__main__':
  args = parse_args()
  preprocess_pdbbind(args.pdbbind_dir)
