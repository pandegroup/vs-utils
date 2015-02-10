"""
Tests for molecule_net.
"""
import os
import unittest

from pande_gas.utils.molecule_net import PcbaJsonParser, PcbaXmlParser


class PcbaParserBase(object):
  """
  Base class for PcbaJsonParser and PcbaXmlParser tests.
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
    self.multiple_target = PcbaJsonParser(
      os.path.join(self.data_dir, 'data/aid429.json'))
    self.gzip_parser = PcbaJsonParser(
      os.path.join(self.data_dir, 'data/aid490.json.gz'))
    self.target_keys = ['name', 'mol_id', 'molecule_type', 'organism']


class TestPcbaXmlParser(unittest.TestCase, PcbaParserBase):
  """
  Tests for PcbaXmlParser.
  """
  def setUp(self):
    """
    Set up tests.
    """
    self.data_dir = os.path.split(os.path.realpath(__file__))[0]
    self.target_keys = ['name', 'mol_id', 'molecule_type']
    self.parser = PcbaXmlParser(
      os.path.join(self.data_dir, 'data/aid490.xml'))
    self.no_target = PcbaXmlParser(
      os.path.join(self.data_dir, 'data/aid1.xml'))
    self.multiple_target = PcbaXmlParser(
      os.path.join(self.data_dir, 'data/aid429.xml'))
    self.gzip_parser = PcbaXmlParser(
      os.path.join(self.data_dir, 'data/aid490.xml.gz'))
