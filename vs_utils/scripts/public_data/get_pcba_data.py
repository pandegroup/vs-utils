"""
Extract PCBA data.

One dataframe is created for each file. The AID and target are associated with
each data point, so routing can be done on a per-point basis using either
field.

Configuration File
------------------
The configuration file is a CSV file whose column headers correspond to columns
that will appear in the saved dataframes. The file must contain an 'aid' column
listing the PCBA assay IDs (AIDs) from which data will be extracted.

Here's an example configuration file:

> aid,target,potency,hill_slope
> 998,757912,Potency,Fit_HillSlope

Running the script with this configuration file will generate a output file
'aid998-gi757912-data.pkl.gz' containing a dataframe with columns ['aid',
'target', 'potency', 'hill_slope', 'efficacy', 'phenotype', 'sid', 'outcome'].

The 'potency' and 'hill_slope' columns will be populated from the 'Potency' and
'Fit_Hillslope' columns in the original data, respectively. The 'aid' and
'target' fields do not match columns in the assay data, so they are considered
constants and will be the same for each row of the dataframe.

Columns are added for fields that are standard for PCBA data, such as a column
to track SIDs ('sid') and categorical activity outcomes ('outcome').
Additionally, columns are added when commonly-occurring fields are recognized
(to simplify writing the configuration file). In this example, 'phenotype' and
'efficacy' columns are added to track the commonly-occurring 'Phenotype' and
'Efficacy' fields.
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
  parser.add_argument('--phenotype', action='store_true',
                      help='Require compound-level phenotype data.')
  parser.add_argument('-s', '--summary',
                      help='Filename for summary information.')
  return parser.parse_args(input_args)


def main(dirs, config_filename, summary_filename=None, with_aid=True,
         with_target=True, phenotype=False):
  aids = set()
  targets = set()
  total = 0
  config = pd.read_csv(config_filename)
  summary = []
  if 'aid' not in config.columns:
    raise ValueError('Configuration file must contain "aid" column.')
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
      try:
        parser = PcbaDataExtractor(filename, this_config, with_aid=with_aid)
      except NotImplementedError as e:
        warnings.warn(e.message)
        continue
      if phenotype and 'phenotype' not in parser.config:
        warnings.warn('{} has no phenotype'.format(aid))
        continue
      assert aid == parser.parser.get_aid()  # sanity check for AID match
      aids.add(aid)
      target = parser.config.get('target')
      targets.add(target)
      data = parser.get_data()
      total += len(data)

      # save dataframe
      output_filename = 'aid{}-{}-data.pkl.gz'.format(aid, target)
      print '{}\t{}\t{}\t{}'.format(aid, target, output_filename, len(data))
      summary.append({'aid': aid, 'target': target,
                      'filename': output_filename, 'size': len(data)})
      write_pickle(data, output_filename)

  # make sure we found everything
  missing = set(config['aid']).difference(aids)
  if len(missing):
    warnings.warn('Missed AIDs {}'.format(missing))

  # save a summary
  summary = pd.DataFrame(summary)
  if summary_filename is not None:
    write_pickle(summary, summary_filename)
  warnings.warn('Found {} assays for {} targets ({} total data points)'.format(
    len(aids), len(targets), total))

if __name__ == '__main__':
  args = parse_args()
  main(args.dirs, args.config, args.summary, args.with_aid, args.with_target,
       args.phenotype)
