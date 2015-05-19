"""
Utilities to prep molecules for nnscore featurization.
"""

import argparse
import openbabel
import os

def force_partial_charge_computation(mol):
  """Force computation of partial charges for molecule.

  This function uses GetPartialCharge to force computation of the Gasteiger
  partial charges. This is an unfortunate hack, since it looks like the
  python openbabel API doesn't expose the OBGastChrg object which actually
  computes partial charges.

  Parameters
  ----------
  mol: OBMol
    Molecule on which we compute partial charges.
  """
  for obatom in openbabel.OBMolAtomIter(mol):
    obatom.GetPartialCharge()


def hydrogenate_and_compute_partial_charges(input_file, input_format, output_directory):
  """Outputs a hydrogenated pdb and a pdbqt with partial charges.

  Takes an input file in specified format. Generates two outputs:

  -) A pdb file that contains a hydrogenated (at pH 7.4) version of
     original compound.
  -) A pdbqt file that has computed Gasteiger partial charges. This pdbqt
     file is build from the hydrogenated pdb.

  Parameters
  ----------
  input_file: String
    Path to input file.
  input_format: String
    Name of input format.
  output_directory: String
    Path to desired output directory.
  """
  basename = os.path.basename(input_file).split(".")[0]

  hyd_output = os.path.join(output_directory, basename + "_hyd.pdb")
  pdbqt_output = os.path.join(output_directory, basename + "_hyd.pdbqt")

  # Create pdb with hydrogens added
  hydConversion = openbabel.OBConversion()
  hydConversion.SetInAndOutFormats(input_format, "pdb")
  mol = openbabel.OBMol()
  hydConversion.ReadFile(mol, input_file)  
  # AddHydrogens(polaronly, correctForPH, pH)
  mol.AddHydrogens(True, True, 7.4)
  hydConversion.WriteFile(mol, hyd_output)

  # Create a pdbqt file from the hydrogenated pdb above.
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

