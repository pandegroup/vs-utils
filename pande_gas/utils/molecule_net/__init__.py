"""
Utilities for MoleculeNet.
"""
import json
import re
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
        self.tree = json.load(filename)

        # should just be one record per file
        assert len(self.tree['PC_AssayContainer']) == 1

        # move in to the assay description
        self.root = self.tree['PC_AssayContainer'][0]['assay']['descr']

    def get_name(self):
        """
        Get assay name.
        """
        return self.root['name']

    def get_description(self):
        """
        Get assay description.
        """
        return '\n'.join(self.root['description'])

    def get_protocol(self):
        """
        Get assay protocol.
        """
        return '\n'.join(self.root['protocol'])

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
        return self.root['target']


class PcbaXmlParser(object):
    """
    Parser for PubChem BioAssay XML.

    Parameters
    ----------
    filename : str
        Filename.
    """
    def __init__(self, filename):
        self.tree = et.parse(filename)
        self.root = self.tree.getroot()

        # default prefix for all tags
        self.prefix = '{http://www.ncbi.nlm.nih.gov}'

    def find(self, tag, root=None):
        """
        Return a list of the elements with a given tag. Note that this only
        searches the direct children of root.

        Parameters
        ----------
        tag : str
            XML tag.
        root : bool, optional (default False)
            Root of XML tree.
        """
        if root is None:
            root = self.root
        return root.findall(self.prefix + tag)

    def join_children(self, elem):
        """
        Join the text for the children of an element.

        Parameters
        ----------
        elem : Element
            Element.
        """
        text = ''
        for child in elem.getchildren():
            if child.text is not None:
                text += child.text + child.tail
        return text

    def get_name(self):
        """
        Get assay name.
        """
        elem = self.find('PC-AssayDescription_name')
        assert len(elem) == 1
        return elem[0].text

    def get_description(self):
        """
        Get assay description.
        """
        elem = self.find('PC-AssayDescription_description')
        assert len(elem) == 1
        return self.join_children(elem[0])

    def get_protocol(self):
        """
        Get assay protocol.
        """
        elem = self.find('PC-AssayDescription_protocol')
        assert len(elem) == 1
        return self.join_children(elem[0])

    def get_target(self):
        """
        Get assay target.

        Returns
        -------
        target : dict
            A dictionary containing keys for target information types, such
            as 'name', 'mol-id', and 'molecule-type'.
        """
        elem = self.find('PC-AssayDescription_target')
        assert len(elem) == 1
        info = self.find('PC-AssayTargetInfo', elem[0])
        assert len(info) == 1
        target = {}
        for e in info[0].getchildren():
            if not e.text.strip():
                continue  # skip blank entries
            m = re.search('PC-AssayTargetInfo_(.*)', e.tag)
            key = m.groups()[0]
            target[key] = e.text
        return target
