"""
Utilities for MoleculeNet.
"""
import gzip
import json
import os
import re
import warnings
import xml.etree.cElementTree as et


class PcbaJsonParser(object):
  """
  Parser for PubChemBioAssay JSON.

  Parameters
  ----------
  filename : str
      Filename.
  """
  def __init__(self, filename):
    if filename.endswith('.gz'):
      with gzip.open(filename) as f:
        self.tree = json.load(f)
    else:
      with open(filename) as f:
        self.tree = json.load(f)

    # move in to the assay description
    self.root = self.tree['PC_AssaySubmit']['assay']['descr']

  def get_aid(self):
    """
    Get assay ID.
    """
    return int(self.root['aid']['id'])

  def get_name(self):
    """
    Get assay name.
    """
    return self.root['name']

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


class PcbaXmlParser(object):
  """
  Parser for PubChem BioAssay XML.

  Parameters
  ----------
  filename : str
      Filename.
  """
  def __init__(self, filename):
    if filename.endswith('.gz'):
      with gzip.open(filename) as f:
        self.tree = et.parse(f)
    else:
      with open(filename) as f:
        self.tree = et.parse(f)

    # default prefix for all tags
    self.prefix = '{http://www.ncbi.nlm.nih.gov}'

    # get pointers to description and data
    descriptions = self.tree.getroot().iter(
      self.prefix + 'PC-AssayDescription')
    descriptions = list(descriptions)
    assert len(descriptions) == 1  # check for only one assay
    self.desc_root = descriptions[0]
    self.data_root = self.tree.getroot().iter(self.prefix + 'PC-AssayResults')

  def find(self, path, root=None):
    """
    Return a list of the elements with a given tag.

    Parameters
    ----------
    path : list
        XML path. Each element of the list is an XML node to which self.prefix
        will be prepended.
    root : bool, optional (default False)
        Root of XML tree.
    """
    if root is None:
      root = self.desc_root
    path = os.path.join(*[self.prefix + tag for tag in path])
    return root.findall(path)

  def join_children(self, elem):
    """
    Join the text for the children of an element.

    Parameters
    ----------
    elem : Element
        Element.
    """
    text = []
    for child in elem.getchildren():
      if child.text is not None:
        text.append(child.text.strip())
      else:
        text.append('')
    return '\n'.join(text)

  def get_aid(self):
    """
    Get assay ID.
    """
    elem = self.find(['PC-AssayDescription_aid', 'PC-ID', 'PC-ID_id'])
    assert len(elem) == 1
    return int(elem[0].text)

  def get_name(self):
    """
    Get assay name.
    """
    elem = self.find(['PC-AssayDescription_name'])
    assert len(elem) == 1
    return elem[0].text

  def get_description(self):
    """
    Get assay description.
    """
    elem = self.find(['PC-AssayDescription_description'])
    assert len(elem) == 1
    return self.join_children(elem[0])

  def get_protocol(self):
    """
    Get assay protocol.
    """
    elem = self.find(['PC-AssayDescription_protocol'])
    assert len(elem) == 1
    return self.join_children(elem[0])

  def get_target(self):
    """
    Get assay target.

    NOTE: Does not return organism information (requires more complicated
        parsing).

    Returns
    -------
    target : dict
        A dictionary containing keys for target information types, such
        as 'name', 'mol-id', and 'molecule-type'.
    """
    # organism information requires more complicated parsing (use JSON)
    warnings.warn('Does not return organism information.')
    elems = self.find(['PC-AssayDescription_target', 'PC-AssayTargetInfo'])
    if len(elems) == 0:
      return None  # no target
    targets = []
    for elem in elems:
      target = {}
      for e in elem.getchildren():
        if not e.text.strip():
          continue  # skip blank entries
        m = re.search('PC-AssayTargetInfo_(.*)', e.tag)
        key = m.groups()[0]
        key = key.replace('-', '_')  # match json
        target[key] = e.text
      targets.append(target)
    assert len(targets)
    return targets
