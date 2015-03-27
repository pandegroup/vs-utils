"""
Parse config file (see get_pcba_data.py) and assemble meta-information for each
target, such as associated AIDs and PDBs.
"""
import argparse
import pandas as pd

from vs_utils.utils import write_pickle


def parse_args(input_args=None):
  """
  Parse command-line arguments.

  Parameters
  ----------
  input_args : list, optional
    Input arguments. If not provided, defaults to sys.argv[1:].
  """
  parser = argparse.ArgumentParser()
  parser.add_argument('config',
                      help='Configuration file.')
  parser.add_argument('-p', '--pdb',
                      help='Target->PDB associations.')
  parser.add_argument('output',
                      help='Output filename.')
  return parser.parse_args(input_args)


def main(config_filename, output_filename, pdb_filename=None):
  """
  Meta-information consists of a row for each target, with a column rach for
  associated AIDs and PDBs (lists).
  """
  # read target->PDB associations
  pdb = {}
  if pdb_filename is not None:
    with open(pdb_filename) as f:
      for line in f:
        target, code = line.split()
        pdb[target] = code.split(',')  # multiple PDBs can be separated by ','
  config = pd.read_csv(config_filename)

  # get AIDs for each target
  targets = {}
  for _, row in config.iterrows():
    target = row['target']
    try:
      int(target)
      target = 'gi_{}'.format(target)  # add 'gi_' to integer targets
    except ValueError:
      pass
    if target not in targets:
      targets[target] = []
    targets[target].append(row['aid'])

  # construct dataframe
  points = []
  for target, aids in targets.iteritems():
    points.append({'target': target, 'aids': aids, 'pdbs': pdb.get(target)})
  df = pd.DataFrame(points)
  write_pickle(df, output_filename)


if __name__ == '__main__':
  args = parse_args()
  main(args.config, args.output, args.pdb)
