"""
Extract PCBA data.

One dataframe is created for each file. The AID and target are associated with
each data point, so routing can be done on a per-point basis using either
field.
"""
import argparse
import glob
import os
import pandas as pd
import warnings

from vs_utils.utils import write_pickle
from vs_utils.utils.public_data import PcbaDataExtractor

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2015, Stanford University"
__license__ = "BSD 3-clause"


def parse_args(input_args=None):
  """
  Parse command-line arguments.

  Parameters
  ----------
  input_args : list, optional
    Input arguments. If not provided, defaults to sys.argv[1:].
  """
  parser = argparse.ArgumentParser()
  parser.add_argument('dirs', nargs='+',
                      help='Directories containing PCBA JSON files.')
  parser.add_argument('-c', '--config', required=1,
                      help='Configuration file containing assay annotations.')
  parser.add_argument('--no-aid', action='store_false', dest='with_aid',
                      help='Do not include AID with each data point.')
  parser.add_argument('--no-target', action='store_false', dest='with_target',
                      help='Do not include target with each data point.')
  return parser.parse_args(input_args)


def main(dirs, config_filename, with_aid, with_target):
  aids = set()
  targets = set()
  total = 0
  config = pd.read_csv(config_filename)
  assert len(config) == len(pd.unique(config['aid']))
  for this_dir in dirs:
    for filename in glob.glob(os.path.join(this_dir, '*.json.gz')):

      # get AID from filename so we only have to load relevant assays
      aid = int(os.path.basename(filename).split('.')[0])
      if aid not in config['aid'].values:
        continue

      # get configuration for this AID
      this_config = config[config['aid'] == aid].iloc[0]
      if not with_aid and 'aid' in this_config:
        del this_config['aid']
      if not with_target and 'target' in this_config:
        del this_config['target']

      # get data
      parser = PcbaDataExtractor(filename, this_config, with_aid=with_aid)
      assert aid == parser.parser.get_aid()  # sanity check for AID match
      aids.add(aid)
      target = parser.config.get('target')
      targets.add(target)
      data = parser.get_data()
      total += len(data)

      # save dataframe
      output_filename = 'aid{}-{}-data.pkl.gz'.format(aid, target)
      print '{}\t{}\t{}\t{}'.format(aid, target, output_filename, len(data))
      write_pickle(data, output_filename)

  # make sure we found everything
  missing = set(config['aid']).difference(aids)
  if len(missing):
    warnings.warn('Missed AIDs {}'.format(missing))

  # print a summary
  print 'Found {} assays for {} targets ({} total data points)'.format(
    len(aids), len(targets), total)

if __name__ == '__main__':
  args = parse_args()
  main(args.dirs, args.config, args.with_aid, args.with_target)
