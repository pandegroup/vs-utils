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
        Return an iterator over elements with a given tag.

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
        elem = list(self.find('PC-AssayDescription_name'))
        assert len(elem) == 1
        return elem[0].text
