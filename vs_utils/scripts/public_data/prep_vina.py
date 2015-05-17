"""
Prep ligand and receptor molecules for loading into vina and nnscore.
"""

import argparse
import openbabel
import os

def parse_args(input_args=None):
  """Parse command-line arguments."""
  parser = argparse.ArgumentParser()
  parser.add_argument('--input_file', required=1,
                      help='Input PDB File.')
  parser.add_argument('--output_directory', required=1,
                      help='Output Directory.')
  return parser.parse_args(input_args)


def force_partial_charge_computation(mol):
  """Force computation of partial charges for molecule.

  This function uses GetPartialCharge to force computation of the Gasteiger
  partial charges. This is an unfortunate hack, since it looks like the
  python openbabel API doesn't expose the OBGastChrg object which actually
  computes partial charges.

  mol: OBMol
    Molecule on which we compute partial charges.
  """
  for obatom in openbabel.OBMolAtomIter(mol):
    obatom.GetPartialCharge()


def main(input_file, output_directory):
  basename = os.path.basename(input_file).split(".")[0]

  hyd_output = os.path.join(output_directory, basename + "_hyd.pdb")
  pdbqt_output = os.path.join(output_directory, basename + "_hyd.pdbqt")

  # Create pdb with hydrogens added
  hydConversion = openbabel.OBConversion()
  hydConversion.SetInAndOutFormats("pdb", "pdb")
  mol = openbabel.OBMol()
  hydConversion.ReadFile(mol, input_file)  
  # AddHydrogens(polaronly, correctForPH, pH)
  mol.AddHydrogens(True, True, 7.4)
  hydConversion.WriteFile(mol, hyd_output)

  # Create a pdbqt file
  chargeConversion = openbabel.OBConversion()
  chargeConversion.SetInAndOutFormats("pdb", "pdbqt")
  # Make protein rigid
  chargeConversion.AddOption("c", chargeConversion.OUTOPTIONS)
  chargeConversion.AddOption("r", chargeConversion.OUTOPTIONS)
  # Preserve hydrogens
  chargeConversion.AddOption("h", chargeConversion.OUTOPTIONS)
  # Preserve atom indices
  chargeConversion.AddOption("p", chargeConversion.OUTOPTIONS)
  # Preserve atom indices
  chargeConversion.AddOption("n", chargeConversion.OUTOPTIONS)

  mol = openbabel.OBMol()
  chargeConversion.ReadFile(mol, hyd_output)
  force_partial_charge_computation(mol)
  chargeConversion.WriteFile(mol, pdbqt_output)


if __name__ == '__main__':
  args = parse_args()
  main(args.input_file, args.output_directory)
