"""
Utilities for MoleculeNet.
"""
import gzip
import json
import re
import warnings
import pandas as pd


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

    # should just be one record per file
    assert len(self.tree['PC_AssayContainer']) == 1

    # move in to the assay description
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
      if (method == "confirmatory"
          and "counterscreen" in self.get_name().lower()):
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

    def add_dataset(self, filename):
      """
      Adds dataset to internal dataframe.
      """
      parser = PcbaJsonParser(filename)
      data = parser.root
      row = {}
      row["name"] = parser.get_name()
      row["aid"] = parser.get_aid()
      row["activity_outcome_method"] = parser.get_activity_outcome_method()
      row["description"] = parser.get_description()
      row["comment"] = parser.get_comment()
      row["results"] = parser.get_results()
      row["revision"] = parser.get_revision()
      self.df.loc[self.index] = pd.Series(row)

    def to_csv(self, out):
      """
      Writes internal dataframe to provided location as csv.
      """
      self.df.to_csv(out)
