"""
Utilities for public_data.
"""
import gzip
import json
import pandas as pd
import warnings


class PcbaJsonParser(object):
  """
  Parser for PubChemBioAssay JSON.

  Parameters
  ----------
  filename : str
      Filename.
  """
  def __init__(self, filename):
    if filename.endswith(".gz"):
      with gzip.open(filename) as f:
        self.tree = json.load(f)
    elif filename.endswith(".json"):
      with open(filename) as f:
        self.tree = json.load(f)
    else:
      raise ValueError("filename must be of type .json or .json.gz!")

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
      name = field['name']
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
    return df

  def get_selected_data(self, column_mapping, include_aid=False, target=None,
                        phenotype=None):
    """
    Get a subset of the assay data.

    Parameters
    ----------
    column_mapping : dict
      Mapping between output dataframe column names and assay data column names
      from which to extract data (e.g. 'potency': 'EC50'). Values that are not
      found in the assay data columns will be treated as constants.
    include_aid : bool, optional (default False)
      Annotate each data point with the AID for this assay.
    target : str, optional
      Assay target to include for each data point.
    phenotype : str, optional
      Default phenotype for non-inactive data points (e.g., 'inhibitor').

    Returns
    -------
    A pandas dataframe containing the selected assay data.
    """
    # add standard PCBA columns
    column_mapping['sid'] = 'sid'
    column_mapping['outcome'] = 'outcome'

    # add optional columns
    if include_aid:
      column_mapping['aid'] = self.get_aid()
    if target is not None:
      column_mapping['target'] = target

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
