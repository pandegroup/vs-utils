"""Build data frames for non-PCBA classification datasets."""
import argparse
import pandas as pd

from vs_utils import utils

from vs_utils.utils.rdkit_utils import serial


def parse_args(input_args=None):
  """Parse command-line arguments."""
  parser = argparse.ArgumentParser()
  parser.add_argument('--assay', required=1,
                      help='Assay ID.')
  parser.add_argument('--target', required=1,
                      help='Assay target.')
  parser.add_argument('-a', '--actives', required=1,
                      help='File containing actives.')
  parser.add_argument('-d', '--decoys', required=1,
                      help='File containing decoys.')
  parser.add_argument('--no-assay', action='store_false', dest='with_assay',
                      help='Do not include AID with each data point.')
  parser.add_argument('--no-target', action='store_false', dest='with_target',
                      help='Do not include target with each data point.')
  parser.add_argument('-j', '--join', action='store_true',
                      help='Join duplicated molecules.')
  parser.add_argument('--phenotype',
                      help='Phenotype for actives in this assay.')
  parser.add_argument('-o', '--output',
                      help='Output filename.')
  parser.add_argument('--mols',
                      help='Filename to write unique molecules.')
  return parser.parse_args(input_args)


def get_rows(reader, outcome, assay, target, with_aid=True, with_target=True,
             phenotype=None, join=False):
  """Get a row for each data point."""
  rows = []
  cids = set()
  mols = set()
  for mol in reader:
    row = {'id': mol.GetProp('_Name'), 'outcome': outcome}
    if join and row['id'] in cids:
      continue
    cids.add(row['id'])
    mols.add(mol)
    if with_aid:
      row['assay'] = assay
    if with_target:
      row['target'] = target
    if phenotype is not None:
      row['phenotype'] = phenotype
    rows.append(row)
  return rows, mols


def main(active_filename, decoy_filename, assay, target, with_assay=True,
         with_target=True, phenotype=None, join=False, output_filename=None,
         unique_filename=None):
  rows = []
  mols = []
  for outcome, filename in zip(['active', 'inactive'],
                               [active_filename, decoy_filename]):
    this_phenotype = phenotype
    if outcome == 'inactive' and phenotype is not None:
      this_phenotype = 'inactive'
    with serial.MolReader().open(filename) as reader:
      this_rows, this_mols = get_rows(reader, outcome, assay, target,
                                      with_assay, with_target, this_phenotype,
                                      join)
      rows.extend(this_rows)
      mols.extend(this_mols)
  assert len(rows) == len(mols)

  df = pd.DataFrame(rows)
  if output_filename is None:
    output_filename = '{}-{}-data.pkl.gz'.format(aid, target)
  print '{}\t{}\t{}\t{}'.format(assay, target, output_filename, len(df))
  utils.write_pickle(df, output_filename)

  if unique_filename is not None:
    with serial.MolWriter().open(unique_filename) as writer:
      writer.write(mols)

if __name__ == '__main__':
  args = parse_args()
  main(args.actives, args.decoys, args.assay, args.target, args.with_assay,
       args.with_target, args.phenotype, args.join, args.output, args.mols)
