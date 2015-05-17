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
  parser.add_argument('--add-hydrogens', action='store_true',
                      dest='add_hydrogens', help='Add hydrogens at pH 7.4')
  parser.add_argument('--gen-pdbqt', action='store_true', dest='gen_pdbqt',
                      help='Generate pdbqt from input.')
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
    print obatom.GetPartialCharge()


def main(input_file, add_hydrogens, gen_pdbqt):
  hydConversion = openbabel.OBConversion()
  if add_hydrogens:
    # Create pdb with hydrogens added
    hydConversion.SetInAndOutFormats("pdb", "pdb")

    mol = openbabel.OBMol()
    hydConversion.ReadFile(mol, input_file)  
    print "Molecule has %d atoms." % mol.NumAtoms()
    # AddHydrogens(polaronly, correctForPH, pH)
    mol.AddHydrogens(True, True, 7.4)
    print "Molecule has %d atoms." % mol.NumAtoms()
    hydConversion.WriteFile(mol, "prgr_hyd.pdb")

  if gen_pdbqt:
    # Create a pdbqt file
    chargeConversion = openbabel.OBConversion()
    chargeConversion.SetInAndOutFormats("pdb", "pdbqt")
    chargeConversion.AddOption("c", chargeConversion.OUTOPTIONS)
    chargeConversion.AddOption("r", chargeConversion.OUTOPTIONS)
    chargeConversion.AddOption("h", chargeConversion.OUTOPTIONS)
    chargeConversion.AddOption("p", chargeConversion.OUTOPTIONS)
    chargeConversion.AddOption("n", chargeConversion.OUTOPTIONS)

    mol = openbabel.OBMol()
    chargeConversion.ReadFile(mol, "prgr_hyd.pdb")
    chargeConversion.WriteFile(mol, "prgr_hyd.pdbqt")


if __name__ == '__main__':
  args = parse_args()
  main(args.input_file, args.add_hydrogens, args.gen_pdbqt)
