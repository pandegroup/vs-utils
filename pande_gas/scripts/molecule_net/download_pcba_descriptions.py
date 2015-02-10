#!/usr/bin/env python
"""
Download PCBA assay descriptions.
"""
import os
import argparse
import gzip
import time

from pubchem_utils import PubChem


def parse_args(input_args=None):
  """
  Parse command-line arguments.

  Parameters
  ----------
  input_args : list, optional
    Input arguments. If not provided, defaults to sys.argv[1:].
  """
  parser = argparse.ArgumentParser()
  parser.add_argument("input",
                      help="Input file containing AIDs.")
  parser.add_argument("-f", "--format", choices=['json', 'xml'],
                      required=True)
  parser.add_argument("--out",
                      help="Output directory for generated files.",
                      required=True)
  parser.add_argument("-n", "--num_files",
                      default=10, type=int,
                      help="Max number of files to download.")
  return parser.parse_args(input_args)


def main(filename, output_format, output_dir, max_num_files):
  """
  Download PCBA JSON descriptions.

  Parameters
  ----------
  filename : str
    Filename containing AIDs.
  output_format : str (default 'json')
    Output file format.
  output_dir : str
    Output directory where downloaded files will be written.
  max_num_files: int
    Maximum Number of file
  """
  if not os.path.isdir(output_dir):
    raise ValueError("%s is not a valid output directory!" % output_dir)
  if filename.endswith(".gz"):
    f = gzip.open(filename)
  elif filename.endswith(".txt"):
    f = open(filename)
  try:
    engine = PubChem()
    count = 0
    start = time.time()
    for line in f:
      if "AID" not in line:
        continue
      # The line is of form 'AID: XXXX' where XXXX is the AID id.
      print line
      aid = int(line.split(":")[1])
      out_file = os.path.join(
          output_dir, 'aid{}.{}'.format(aid, output_format))
      data = engine.get_assay_description(aid, output_format)
      with open(out_file, 'wb') as f:
          f.write(data)
      count += 1
      if count >= max_num_files:
        break
  finally:
    f.close()
    end = time.time()
    print "Elapsed Time: " + str(end-start)

if __name__ == '__main__':
  args = parse_args()
  main(args.input, args.format, args.out, args.num_files)
