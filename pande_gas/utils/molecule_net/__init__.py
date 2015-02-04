"""
Utilities for MoleculeNet.
"""
import xml.etree.cElementTree as et


class AssayXMLParser(object):
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

    def find(self, tag):
        """
        Return a list of the elements with a given tag.

        Parameters
        ----------
        tag : str
            XML tag.
        """
        return self.root.findall(self.prefix + tag)

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
        description = ''
        for child in elem[0].getchildren():
            if child.text is not None:
                description += child.text + child.tail
        return description

    def get_protocol(self):
        """
        Get assay protocol.
        """

    def get_target(self):
        """
        Get assay target.
        """
