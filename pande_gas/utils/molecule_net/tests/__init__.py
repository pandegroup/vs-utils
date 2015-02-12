"""
Tests for molecule_net.
"""
import os
import unittest

from pande_gas.utils.molecule_net import PcbaJsonParser, PcbaPandasHandler


class PcbaParserBase(object):
  """
  Base class for PcbaJsonParser tests.
  """
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

  # TODO(rbharath): add a test for get_results. The problem is that get_results
  # is a dict with many fields. Which if any of these do we want to recognize
  # explicitly?

  def test_get_revision(self):
    """
    Test parsing of revision.
    """
    revision = 1
    assert self.no_target.get_revision() == revision


class TestPcbaJsonParser(unittest.TestCase, PcbaParserBase):
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
    self.target_keys = ['name', 'mol_id', 'molecule_type', 'organism']


class TestPcbaPandasHandler(unittest.TestCase):
  """
  Tests for PcbaPandasHandler.
  """
  def setUp(self):
    self.handler = PcbaPandasHandler()
    self.data_dir = os.path.split(os.path.realpath(__file__))[0]

  def test_add_dataset(self):
    """
    Test for adding datasets.
    """
    self.handler.add_dataset(
        os.path.join(self.data_dir, "data/aid1.json"))
    num_rows = len(self.handler.df.index)
    assert num_rows == 1

  def test_to_csv(self):
    """
    Test for witing to csv
    """
    self.handler.add_dataset(
        os.path.join(self.data_dir, "data/aid1.json"))
    self.handler.to_csv("/usr/local/google/home/bramsundar/pande-gas/pande_gas/out.csv")
