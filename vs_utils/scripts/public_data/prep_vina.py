"""
Prep ligand and receptor molecules for loading into vina and nnscore.
"""

import argparse
import openbabel
import os
from vs_utils.scripts.public_data.nnscore_utilities import hydrogenate_and_compute_partial_charges

def parse_args(input_args=None):
  """Parse command-line arguments."""
  parser = argparse.ArgumentParser()
  parser.add_argument('--input-file', required=1,
                      help='Input File.')
  parser.add_argument('--output-directory', required=1,
                      help='Output Directory.')
  parser.add_argument('--input-format', required=1,
                      help='File format of input file.')
  return parser.parse_args(input_args)

def main(input_file, input_format, output_directory):
  hydrogenate_and_compute_partial_charges(input_file, input_format,
      output_directory)


if __name__ == '__main__':
  args = parse_args()
  main(args.input_file, args.input_format, args.output_directory)
