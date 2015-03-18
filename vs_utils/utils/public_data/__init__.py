"""
Utilities for public_data.
"""
import gzip
try:
  import ujson as json
except ImportError:
  import json
import numpy as np
import pandas as pd
import warnings


def read_json(filename):
  """
  Read a JSON file.

  Parameters
  ----------
  filename : str
    Filename. Must be of type .json or .json.gz.
  """
  if filename.endswith('json.gz'):
    with gzip.open(filename) as f:
      tree = json.load(f)
  elif filename.endswith('.json'):
    with open(filename) as f:
      tree = json.load(f)
  else:
    raise ValueError('Filename must be of type .json or .json.gz.')
  return tree


def read_sid_cid_map(filename):
  """
  Read SID->CID map.

  Parameters
  ----------
  filename : str
    SID->CID map.
  """
  if filename.endswith('.gz'):
    f = gzip.open(filename)
  else:
    f = open(filename)
  try:
    sid_cid = {}
    for line in f:
      sid, cid = line.split()
      assert int(sid) not in sid_cid
      sid_cid[int(sid)] = int(cid)
    return sid_cid
  finally:
    f.close()


class PcbaJsonParser(object):
  """
  Parser for PubChemBioAssay JSON.

  Parameters
  ----------
  filename : str
      Filename.
  """
  def __init__(self, filename):
    self.tree = read_json(filename)
    self.data = None

    # move in to the assay description
    try:
        # FTP format
        self.root = self.tree['PC_AssaySubmit']['assay']['descr']
    except KeyError:
        # REST format
        # should just be one record per file
        assert len(self.tree['PC_AssayContainer']) == 1
        self.root = self.tree['PC_AssayContainer'][0]['assay']['descr']

  def get_name(self):
    """
    Get assay name.
    """
    return self.root['name']

  def get_aid(self):
    """
    Get assay AID.
    """
    return self.root["aid"]["id"]

  def get_activity_outcome_method(self):
    """
    Get activity outcome method.
    """
    #
    if "activity_outcome_method" in self.root:
      method = self.root["activity_outcome_method"]
      if "counter" in self.get_name().lower():
        method = "counterscreen"
      return method
    else:
      return None

  def get_description(self):
    """
    Get assay description.
    """
    if isinstance(self.root['description'], list):
      return '\n'.join(
        [line.strip() for line in self.root['description']])
    else:
      return self.root['description']

  def get_protocol(self):
    """
    Get assay protocol.
    """
    if isinstance(self.root['protocol'], list):
      return '\n'.join([line.strip() for line in self.root['protocol']])
    else:
      return self.root['protocol']

  def get_target(self):
    """
    Get assay target.

    TODO: Decide which fields are important. We may be able to match
        targets by mol-id.

    Returns
    -------
    target : dict
        A dictionary containing keys for target information types, such
        as 'name', 'mol-id', and 'molecule-type'.
    """
    if 'target' in self.root:
      return self.root['target']
    else:
      return None

  def get_comment(self):
    """
    Get assay comment.
    """
    if "comment" in self.root:
      if isinstance(self.root["comment"], list):
        return "\n".join([line.strip() for line in self.root["comment"]])
      else:
        return self.root["comment"]
    else:
      return None

  def get_results(self):
    """
    Get Assay result fields.
    """
    if "results" in self.root:
      return self.root["results"]
    else:
      return None

  def get_revision(self):
    """
    Get assay revision.
    """
    if "revision" in self.root:
      return self.root["revision"]
    else:
      return None

  def get_result_names(self, from_tid=False):
    """
    Get column names for assay-specific result fields.

    Parameters
    ----------
    from_tid : bool, optional (default False)
      Return a dict mapping TIDs to field names. If False, returns a list of
      field names.
    """
    tids = {}
    names = []
    for field in self.get_results():
      name = field['name'].strip()  # clean up extra whitespace
      if name in names:
        warnings.warn(
          'Duplicated field "{}" in AID {}'.format(name, self.get_aid()))
      tids[field['tid']] = name
      names.append(name)
    if from_tid:
      return tids
    else:
      return names

  def get_data(self):
    """
    Get assay data in a Pandas dataframe.
    """
    if self.data is not None:
      return self.data
    try:
      data = self.tree['PC_AssaySubmit']['data']
    except KeyError:
      return None

    # populate dataframe
    tids = self.get_result_names(from_tid=True)
    series = []
    for dp in data:
      point = {}
      for key, value in dp.iteritems():
        if key == 'data':  # assay-specific fields
          for col in value:
            col_name = tids[col['tid']]
            assert len(col['value']) == 1
            for col_value in col['value'].itervalues():
              point[col_name] = col_value
        else:  # generic fields
          point[key] = value
      series.append(point)
    df = pd.DataFrame(series)
    assert len(df) == len(data)
    self.data = df
    return df

  def get_selected_data(self, column_mapping, with_aid=False, phenotype=None):
    """
    Get a subset of the assay data.

    Parameters
    ----------
    column_mapping : dict
      Mapping between output dataframe column names and assay data column names
      from which to extract data (e.g. 'potency': 'EC50'). Values that are not
      found in the assay data columns will be treated as constants.
    with_aid : bool, optional (default False)
      Annotate each data point with the AID for this assay.
    phenotype : str, optional
      Default phenotype for non-inactive data points (e.g., 'inhibitor').

    Returns
    -------
    A pandas dataframe containing the selected assay data.
    """
    # make a copy of the column mapping
    column_mapping = column_mapping.copy()

    # add standard PCBA columns
    column_mapping['sid'] = 'sid'
    column_mapping['outcome'] = 'outcome'

    # add optional columns
    if with_aid:
      column_mapping['aid'] = self.get_aid()
    elif 'aid' in column_mapping:
      warnings.warn('column_mapping contains "aid"')

    # get selected columns
    data = self.get_data()
    if data is None:
      return None
    old_cols, new_cols = [], []
    constants = {}
    for new_col, old_col in column_mapping.iteritems():
      if old_col not in data.columns:
        constants[new_col] = old_col
      else:
        new_cols.append(new_col)
        old_cols.append(old_col)
    df = data[old_cols]  # get selected columns from assay data
    df.columns = new_cols  # rename columns to match column_mapping
    for new_col, value in constants.iteritems():
      df.insert(0, new_col, value)  # add constant-valued columns

    # process phenotypes
    if phenotype is not None and 'phenotype' not in column_mapping:
      df.insert(len(df.columns), 'phenotype', phenotype)
      df.loc[df['outcome'] == 'inactive', 'phenotype'] = 'inactive'

    return df


class PcbaPandasHandler(object):
    """
    Writes data from PCBA into pandas dataframes.

    Parameters
    ----------
    """
    def __init__(self):
      self.index = 0
      self.df = pd.DataFrame(
          columns=["name", "aid", "activity_outcome_method",
                   "description", "comment", "results", "revision"])
      self.df['aid'] = self.df['aid'].astype(int)  # force AID to int

    def add_dataset(self, filename):
      """
      Adds dataset to internal dataframe.
      """
      parser = PcbaJsonParser(filename)
      row = {}
      row["name"] = parser.get_name()
      row["aid"] = parser.get_aid()
      row["activity_outcome_method"] = parser.get_activity_outcome_method()
      row["description"] = parser.get_description()
      row["comment"] = parser.get_comment()
      row["results"] = parser.get_results()
      row["revision"] = parser.get_revision()
      self.df.loc[self.index] = pd.Series(row)
      self.index += 1  # increment index

    def get_dataset(self, index):
      """
      Fetches information for a particular dataset by index.
      """
      return self.df.loc[index]

    def to_csv(self, out):
      """
      Writes internal dataframe to provided location as csv.
      """
      self.df.to_csv(out)


class PcbaDataExtractor(object):
  """
  Extract selected data from PCBA assay data.

  Parameters
  ----------
  filename : str
    PCBA JSON assay data file.
  config : dict or pd.Series
    Mapping between output dataframe column names and assay data column names
    from which to extract data (e.g. 'potency': 'EC50'). Values that are not
    found in the assay data columns will be treated as constants.
  with_aid : bool, optional (default False)
    Annotate each data point with the AID for this assay.
  """
  def __init__(self, filename, config, with_aid=True):
    self.filename = filename
    self.parser = PcbaJsonParser(filename)
    self.with_aid = with_aid
    self.config = config
    self.phenotype = None  # default phenotype for this assay
    self._check_config()  # check configuration

  def get_data(self, lower=True, sid_cid=None):
    """
    Get selected data from the assay.

    Parameters
    ----------
    lower : bool, optional (default True)
      Lowercase string fields for consistency.
    sid_cid : dict, optional
      SID->CID mapping. If provided, adds a 'cid' column to the dataframe.
    """
    data = self.parser.get_selected_data(self.config, with_aid=self.with_aid,
                                         phenotype=self.phenotype)

    # map SIDs to CIDs
    if sid_cid is not None:
      cids = [sid_cid.get(sid) for sid in data['sid'].values]
      data['cid'] = pd.Series(cids, index=data.index)

    # lowercase string fields for consistency
    if lower:
      for col, dtype in data.dtypes.iteritems():
        if dtype == np.dtype('object'):
          data.loc[:, col] = data[col].str.lower()
    return data

  def _check_config(self):
    """
    Check the column mapping and other configuration parameters.
    """
    # make a copy of the config
    config = self.config.copy()

    # remove null columns
    for key, value in config.iteritems():
      if pd.isnull(value):
        del config[key]

    # check for some common column names
    columns = self.parser.get_result_names()
    if 'Potency' in columns:
      config['potency'] = 'Potency'
    if 'Efficacy' in columns:
      config['efficacy'] = 'Efficacy'
    if 'Phenotype' in columns:
      config['phenotype'] = 'Phenotype'

    # target cleanup
    # add gi: prefix to integer targets
    if 'target' in config:
      try:
        int(config['target'])
        config['target'] = 'gi_{}'.format(config['target'])
      except ValueError:
        pass

    # phenotype handling
    # either a column name or a default annotation
    phenotypes = {'in': 'inhibitor', 'ac': 'activator', 'ag': 'activator',
                  'an': 'inhibitor', 'ia': 'inhibitor', 'pam': 'PAM',
                  'nam': 'NAM', '?': None}
    phenotype = config.get('phenotype')
    if phenotype and phenotype not in columns:
      del config['phenotype']  # remove from column mapping
      if phenotype in phenotypes.iterkeys():
        phenotype = phenotypes[phenotype]  # map to full name
      if phenotype not in phenotypes.itervalues():
        raise NotImplementedError(
          'Unrecognized phenotype "{}"'.format(phenotype))
      self.phenotype = phenotype  # set default phenotype for this assay

    self.config = config
