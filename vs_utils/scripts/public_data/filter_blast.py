#!/usr/bin/env python
"""
Filter BLAST results. Threshold results based on various metrics such as
expectation value, sequence identity, and query or subject coverage.

Some of these metrics are not part of the standard BLAST output and must be
specifically requested (see blast.sh).
"""
import argparse
import numpy as np
import pandas as pd


def parse_args(input_args=None):
  """
  Parse command-line arguments.

  Parameters
  ----------
  input_args : list, optional
    Input arguments. If not provided, defaults to sys.argv[1:].
  """
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input', required=1,
                      help='BLAST results.')
  parser.add_argument('-o', '--output', required=1,
                      help='Output filename.')
  parser.add_argument('--e-value', type=float, default=1,
                      help='Maximum E-value.')
  parser.add_argument('--identity', type=float, default=90,
                      help='Minimum sequence identity (%).')
  parser.add_argument('--query-coverage', type=float, default=90,
                      help='Minimum query coverage (%).')
  parser.add_argument('--subject-coverage', type=float, default=90,
                      help='Minimum subject coverage (%).')
  return parser.parse_args(input_args)


def main(input_filename, output_filename, e_value=1, identity=90,
         query_coverage=90, subject_coverage=90):
  df = pd.read_table(input_filename)
  df['gi'] = pd.DataFrame(df['query'].str.split('|').tolist())[1]
  df['pdb'] = pd.DataFrame(df['subject'].str.split('_').tolist())[0]
  df['query_coverage'] = 100 * np.true_divide(df['alignment_length'],
                                              df['query_length'])
  df['subject_coverage'] = 100 * np.true_divide(df['alignment_length'],
                                                df['subject_length'])
  filtered = df[(df.e_value <= e_value)
                & (df.identity >= identity)
                & (df.query_coverage >= query_coverage)
                & (df.subject_coverage >= subject_coverage)]
  print '{} -> {}'.format(df.shape[0], filtered.shape[0])
  filtered.to_csv(output_filename)

  # get unique target/PDB pairings
  pairings = {}
  for gi, pdb in zip(filtered['gi'], filtered['pdb']):
    if gi not in pairings:
      pairings[gi] = set()
    pairings[gi].add(pdb)
  for gi, pdbs in pairings.iteritems():
    for pdb in pdbs:
      print '{}\t{}'.format(gi, pdb)

if __name__ == '__main__':
  args = parse_args()
  main(args.input, args.output, args.e_value, args.identity,
       args.query_coverage, args.subject_coverage)
