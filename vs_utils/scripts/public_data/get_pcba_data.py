"""
Extract PCBA data.

One dataframe is created for each file. The AID and target are associated with
each data point, so routing can be done on a per-point basis using either
field.
"""
import argparse
import glob
import numpy as np
import os
import pandas as pd

from vs_utils.utils import write_pickle
from vs_utils.utils.public_data import PcbaJsonParser

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
  parser.add_argument('--no-aid', action='store_false', dest='include_aid',
                      help='Do not include AID with each data point.')
  parser.add_argument('--no-target', action='store_false',
                      dest='include_target',
                      help='Do not include target with each data point.')
  parser.add_argument('--config', required=1,
                      help='Configuration file containing assay annotations.')
  return parser.parse_args(input_args)


def main(dirs, config_filename, include_aid, include_target):
  aids = set()
  targets = set()
  total = 0
  phenotypes = {'in': 'inhibitor', 'ac': 'activator', 'ag': 'activator',
                'an': 'inhibitor', 'ia': 'inhibitor'}
  config = pd.read_csv(config_filename)
  assert len(config) == len(pd.unique(config['aid']))
  for this_dir in dirs:
    for filename in glob.glob(os.path.join(this_dir, '*.json.gz')):
      file_aid = int(os.path.basename(filename).split('.')[0])
      if file_aid not in config['aid'].values:
        continue  # only load relevant assays
      parser = PcbaJsonParser(filename)
      aid = parser.get_aid()
      aids.add(aid)
      this_config = config[config['aid'] == aid].iloc[0]

      # check common column names
      columns = parser.get_result_names()
      if 'Potency' in columns:
        this_config['potency'] = 'Potency'
      if 'Efficacy' in columns:
        this_config['efficacy'] = 'Efficacy'
      if 'Phenotype' in columns:
        this_config['phenotype'] = 'Phenotype'

      # target cleanup
      target = None
      if include_target and 'target' in this_config:
        target = this_config['target']
        try:
          int(this_config['target'])  # add gi: prefix to integer targets
          target = 'gi{}'.format(this_config['target'])
          this_config['target'] = target
        except ValueError:
          pass
      targets.add(target)

      # phenotype handling
      # either a column name or a default annotation
      phenotype = this_config.get('phenotype', None)
      if phenotype and phenotype not in parser.get_result_names():
        del this_config['phenotype']  # remove from column mapping
        if phenotype in phenotypes.iterkeys():
          phenotype = phenotypes[phenotype]  # map to full name
        if phenotype not in phenotypes.itervalues():
          raise NotImplementedError(
            'Unrecognized phenotype "{}"'.format(phenotype))
      else:
        phenotype = None

      data = parser.get_selected_data(this_config, include_aid=include_aid,
                                      target=target, phenotype=phenotype)

      # lowercase string fields for consistency
      for col, dtype in data.dtypes.iteritems():
        if dtype == np.dtype('object'):
          data[col] = data[col].str.lower()

      output_filename = 'aid{}-{}-data.pkl.gz'.format(aid, target)
      total += len(data)
      print '{}\t{}\t{}\t{}'.format(aid, target, output_filename, len(data))
      write_pickle(data, output_filename)
  assert len(aids) == len(config)  # make sure we found everything
  print 'Found {} assays for {} targets ({} total data points)'.format(
    len(aids), len(targets), total)

if __name__ == '__main__':
  args = parse_args()
  main(args.dirs, args.config, args.include_aid, args.include_target)
