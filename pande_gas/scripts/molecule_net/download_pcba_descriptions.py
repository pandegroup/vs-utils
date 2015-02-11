#!/usr/bin/env python
"""
Download PCBA assay descriptions.
"""
import argparse
import gzip
import json
import os

from pubchem_utils import PubChem

from pande_gas.utils import write_pickle


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
  parser.add_argument("out",
                      help="Output directory for descriptions.")
  parser.add_argument('-np', '--n_jobs', type=int, default=1,
                      help='Number of parallel jobs.')
  return parser.parse_args(input_args)


def read_aids(filename):
  """
  Read AIDs from file.

  Parameters
  ----------
  filename : str
    Filename containing AIDs.
  """
  if filename.endswith('.gz'):
    f = gzip.open(filename)
  else:
    f = open(filename)
  try:
    aids = []
    for line in f:
      if line.strip():
        aids.append(int(line))
  finally:
    f.close()
  return aids


def main(filename, output_dir, n_jobs=1):
  """
  Download PCBA JSON descriptions.

  Parameters
  ----------
  filename : str
    Filename containing AIDs.
  output_dir : str
    Output directory for assay descriptions.
  n_jobs : int (default 1)
    Number of parallel jobs.
  """
  aids = read_aids(filename)
  print 'Downloading JSON descriptions for {} assays...'.format(len(aids))
  engine = PubChem()
  descriptions = engine.get_assay_descriptions(aids, n_jobs=n_jobs)
  
  # create output directory if it doesn't exist
  try:
    os.mkdir(output_dir)
  except OSError:
    pass

  # write descriptions to individual files
  for description in descriptions:
    aid = description['aid']['id']
    assert aid
    
    # nest description so it looks like PUG REST output
    nested = {'PC_AssayContainer': [{'assay': {'descr': description}}]}

    with gzip.open(os.path.join(output_dir, 
                   'aid{}.json.gz'.format(aid)), 'wb') as f:
      json.dump(nested, f)

if __name__ == '__main__':
  args = parse_args()
  main(args.input, args.out, args.n_jobs)
