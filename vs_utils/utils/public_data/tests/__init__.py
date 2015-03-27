"""
Tests for public_data.
"""
import gzip
import os
import pandas as pd
import tempfile
import unittest
import csv

from vs_utils.utils.public_data import (PcbaDataExtractor, PcbaJsonParser,
                                        PcbaPandasHandler, read_json,
                                        read_sid_cid_map)


def test_read_json():
  """
  Test read_json.
  """
  data = read_json(os.path.join(os.path.split(os.path.realpath(__file__))[0],
                                'data/aid490.json'))
  assert data['PC_AssaySubmit']['assay']['descr']['aid']['id'] == 490
  data = read_json(os.path.join(os.path.split(os.path.realpath(__file__))[0],
                                'data/aid490.json.gz'))
  assert data['PC_AssaySubmit']['assay']['descr']['aid']['id'] == 490


def test_read_sid_cid_map():
  """
  Test read_sid_cid_map.
  """
  f, filename = tempfile.mkstemp(suffix='.txt')
  os.close(f)
  g, gilename = tempfile.mkstemp(suffix='.txt.gz')
  os.close(g)
  try:
    with open(filename, 'wb') as f:
      f.write('123456\t7890\n')
    with gzip.open(gilename, 'wb') as g:
      g.write('123456\t7890\n')
    assert read_sid_cid_map(filename) == {123456: 7890}
    assert read_sid_cid_map(gilename) == {123456: 7890}
  finally:
    os.remove(filename)
    os.remove(gilename)


class TestPcbaJsonParser(unittest.TestCase):
  """
  Tests for PcbaJsonParser.
  """
  def setUp(self):
    """
    Set up tests.
    """
    self.data_dir = os.path.split(os.path.realpath(__file__))[0]
    self.parser = PcbaJsonParser(
      os.path.join(self.data_dir, 'data/aid490.json'))
    self.no_target = PcbaJsonParser(
      os.path.join(self.data_dir, 'data/aid1.json'))
    self.confirmatory = self.no_target
    self.multiple_target = PcbaJsonParser(
      os.path.join(self.data_dir, 'data/aid429.json'))
    self.gzip_parser = PcbaJsonParser(
      os.path.join(self.data_dir, 'data/aid490.json.gz'))
    self.rest_parser = PcbaJsonParser(
      os.path.join(self.data_dir, 'data/aid1-rest.json'))
    self.data_parser = PcbaJsonParser(
      os.path.join(self.data_dir, 'data/540325.json.gz'))
    self.target_keys = ['name', 'mol_id', 'molecule_type', 'organism']

  def test_get_aid(self):
    """
    Test get_aid.
    """
    assert self.parser.get_aid() == 490

  def test_get_name(self):
    """
    Test get_name.
    """
    name = ('Literature data for small-molecule inhibitors of ' +
            'Influenza_A_virus_(A_Tokyo_3_67(H2N2))')
    assert self.parser.get_name() == name

  def test_get_description(self):
    """
    Test get_description.
    """
    desc = ('This assay contains in vitro affinity data extracted from ' +
            'the literature for compounds tested against Influenza A ' +
            'virus (strain A/Tokyo/3/67 H2N2).')
    assert self.parser.get_description() == desc

  def test_get_protocol(self):
    """
    Test get_protocol.
    """
    protocol = ('Various protocols were used; consult BindingDB and/or ' +
                'cited papers for details.\n\nA compound is listed as ' +
                '\"Active\" if IC50<100,000 nanomolar or Ki<100,000 ' +
                'nanomolar.\n\nIf multiple measurements are available ' +
                'for a given compound, the compound is listed as ' +
                '\"Active\" if any of the measurements meet the ' +
                'criterion.')
    assert self.parser.get_protocol() == protocol

  def test_get_target(self):
    """
    Test get_target.
    """
    target = self.parser.get_target()
    assert len(target) == 1
    for key in self.target_keys:
      assert key in target[0]

  def test_get_no_target(self):
    """
    Test get_target on an assay with no target.
    """
    assert self.no_target.get_target() is None

  def test_get_multiple_target(self):
    """
    Test get_target on an assay with multiple targets.
    """
    targets = self.multiple_target.get_target()
    assert len(targets) == 2
    for target in targets:
      for key in self.target_keys:
        assert key in target, key

  def test_parse_gzip(self):
    """
    Test parsing gzipped files.
    """
    name = ('Literature data for small-molecule inhibitors of ' +
            'Influenza_A_virus_(A_Tokyo_3_67(H2N2))')
    assert self.gzip_parser.get_name() == name

  def test_get_activity_outcome_method(self):
    """
    Test parsing of activity_outcome_method.
    """
    method = "confirmatory"
    assert self.confirmatory.get_activity_outcome_method() == method

  def test_get_comment(self):
    """
    Test parsing of comment.
    """
    comment = ("These data are a subset of the data "
               + "from the NCI human tumor cell line screen. "
               + "Compounds are identified by the NCI NSC number. "
               + "In the NCI numbering system, NCI-H23 is panel "
               + "number 1, cell number 1\nBasically compounds  "
               + "with LogGI50 (unit M) less than -6 were considered "
               + "as active. Activity score was based on increasing "
               + "values of -LogGI50.")
    assert self.no_target.get_comment() == comment

  def test_get_results(self):
    """
    Test parsing of results.

    TODO(rbharath): get_results returns a dict with many fields. Which if
    any of these do we want to recognize explicitly?
    """
    results = [
      {u'name': u'loggi50', u'transform': u'log', u'tid': 1, u'type':
       u'float', u'unit': u'm', u'description': [u'Log of the GI50 result, '
                                                 u'unit: M.']},
     {u'name': u'loggi50', u'transform': u'log', u'tid': 2, u'type':
      u'float', u'unit': u'ugml', u'description': [u'Log of the GI50 '
                                                   u'result, unit: ug/mL.']},
     {u'name': u'loggi50', u'transform': u'log', u'sunit': u'v', u'tid': 3,
      u'type': u'float', u'description': [u'Log of the GI50 result, unit: '
                                          u'Volumetric.']},
     {u'tid': 4, u'type': u'int', u'name': u'indngi50', u'description':
      [u'Number of tests averaged for the GI50 value for this NSC and cell '
       u'line.']},
     {u'tid': 5, u'transform': u'log', u'type': u'float', u'name':
      u'stddevgi50', u'description': [u'Standard Deviation of the Log10 of '
                                      u'the GI50 result averaged across all '
                                      u'tests for this NSC and cell line.']},
     {u'name': u'logtgi', u'transform': u'log', u'tid': 6, u'type':
      u'float', u'unit': u'm', u'description': [u'Log of the TGI result, '
                                                u'unit: M.']},
     {u'name': u'logtgi', u'transform': u'log', u'tid': 7, u'type':
      u'float', u'unit': u'ugml', u'description': [u'Log of the TGI result, '
                                                   u'unit: ug/mL.']},
     {u'name': u'logtgi', u'transform': u'log', u'sunit': u'v', u'tid': 8,
      u'type': u'float', u'description': [u'Log of the TGI result, unit: '
                                          u'Volumetric.']},
     {u'tid': 9, u'type': u'int', u'name': u'indntgi', u'description':
      [u'Number of tests averaged for the TGI value for this NSC and cell '
       u'line.']},
     {u'tid': 10, u'transform': u'log', u'type': u'float', u'name':
      u'stddevtgi', u'description': [u'Standard Deviation of the Log10 of '
                                     u'the TGI result averaged across all '
                                     u'tests for this NSC and cell line.']}]
    assert self.no_target.get_results() == results

  def test_get_revision(self):
    """
    Test parsing of revision.
    """
    revision = 1
    assert self.no_target.get_revision() == revision

  def test_read_rest_json(self):
    """
    Test parsing of REST-formatted JSON.
    """
    assert self.rest_parser.get_aid() == 1

  def test_get_data(self):
    """
    Test get_data.
    """
    df = self.data_parser.get_data()
    assert df.shape == (110, 6)
    rows = df[df['sid'] == 24831307]
    assert len(rows) == 1
    assert rows.iloc[0]['Activity at 10 uM'] == 99.68

    # check for no data
    assert self.parser.get_data() is None

  def test_get_selected_data(self):
    """
    Test get_selected_data.
    """
    config = {'activity': 'Activity at 10 uM', 'blah': 'constant'}
    data = self.data_parser.get_selected_data(config, with_aid=True)
    assert 'aid' not in config  # check that the method didn't change config
    assert len(data) == 110
    row = data[data['sid'] == 24831307].iloc[0]
    assert row['aid'] == 540325
    assert row['activity'] == 99.68
    assert row['blah'] == 'constant'

    # check for no data
    assert self.parser.get_selected_data(config) is None

  def test_get_result_names(self):
    """
    Test get_result_names.
    """
    results = {1: 'Ki_max', 2: 'Ki', 3: 'Ki_min', 4: 'IC50_max', 5: 'IC50',
               6: 'IC50_min', 7: 'PubMed Citation (PMID)',
               8: 'WILD TYPE sequence'}
    assert set(self.parser.get_result_names()) == set(results.values())
    assert self.parser.get_result_names(from_tid=True) == results


class TestPcbaPandasHandler(unittest.TestCase):
  """
  Tests for PcbaPandasHandler.
  """
  def setUp(self):
    self.handler = PcbaPandasHandler()
    self.data_dir = os.path.split(os.path.realpath(__file__))[0]
    self.parser = PcbaJsonParser(
      os.path.join(self.data_dir, 'data/aid1.json'))

  def test_add_dataset(self):
    """
    Test for adding datasets.
    """
    self.handler.add_dataset(
        os.path.join(self.data_dir, "data/aid1.json"))
    row = self.handler.get_dataset(0)
    num_rows = len(self.handler.df)
    assert num_rows == 1
    assert row.get("comment") == self.parser.get_comment()
    assert row.get("name") == self.parser.get_name()
    assert row.get("aid") == self.parser.get_aid()

  def test_to_csv(self):
    """
    Test for witing to csv
    """
    # Open and close a tempfile without deleting it. (On some systems
    # cannot open tempfiles if they are already open, so need to close
    # to use elsewhere in test)
    f = tempfile.NamedTemporaryFile(delete=False)
    f.close()
    try:
      self.handler.add_dataset(
          os.path.join(self.data_dir, "data/aid1.json"))
      self.handler.to_csv(f.name)
      with open(f.name, "rb") as csvfile:
        reader = csv.DictReader(csvfile)
        #for row in reader:
        row = reader.next()
        assert row["comment"] == self.parser.get_comment()
        assert row["name"] == self.parser.get_name()
        assert int(row["aid"]) == self.parser.get_aid()
    finally:
      # Delete tempfile
      os.remove(f.name)


class TestPcbaDataExtractor(unittest.TestCase):
  """
  Tests for PcbaDataExtractor.
  """
  def setUp(self):
    """
    Set up tests.
    """
    self.data_dir = os.path.split(os.path.realpath(__file__))[0]
    self.aid998 = os.path.join(self.data_dir, 'data/998.json.gz')
    self.aid540325 = os.path.join(self.data_dir, 'data/540325.json.gz')

  def test_get_data(self):
    """
    Test get_data.
    """
    config = {'target': '1296534', 'r2': 'Fit_R2', 'phenotype': 'in'}
    engine = PcbaDataExtractor(self.aid998, config, with_aid=True)

    # check SID->CID mapping
    data = engine.get_data()
    assert 'cid' not in data.columns
    data = engine.get_data(sid_cid={11110959: 1730})
    assert data[data['sid'] == 11110959].iloc[0]['cid'] == 1730
    assert pd.isnull(data[data['sid'] == 11111313].iloc[0]['cid'])

    # check lowercase
    data = engine.get_data(lower=True)
    assert data[data['sid'] == 11110959].iloc[0]['phenotype'] == 'inhibitor'
    assert data[data['sid'] == 11110959].iloc[0]['potency'] == 4.4668
    data = engine.get_data(lower=False)
    assert data[data['sid'] == 11110959].iloc[0]['phenotype'] == 'Inhibitor'
    assert data[data['sid'] == 11110959].iloc[0]['potency'] == 4.4668

  def test_check_config(self):
    """
    Test _check_config.
    """
    config = {'target': '1296534', 'potency': 'Activity at 10 uM',
              'phenotype': 'in'}

    # without common fields
    engine = PcbaDataExtractor(self.aid540325, config, with_aid=True)
    assert engine.config['target'] == 'gi_1296534'
    assert 'phenotype' not in engine.config
    assert engine.phenotype == 'inhibitor'

    # with common fields
    engine = PcbaDataExtractor(self.aid998, config, with_aid=True)
    assert engine.config['target'] == 'gi_1296534'
    assert 'phenotype' in engine.config
    assert engine.config['phenotype'] == 'Phenotype'
    assert engine.config['potency'] == 'Potency'
    assert engine.config['efficacy'] == 'Efficacy'
