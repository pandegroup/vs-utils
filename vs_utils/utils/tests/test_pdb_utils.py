"""
Tests for pdb_utils.
"""
from cStringIO import StringIO
import unittest

from vs_utils.utils import pdb_utils


class TestPDBUtils(unittest.TestCase):
    """
    Tests for pdb_utils.
    """
    def setUp(self):
        """
        Set up tests.
        """
        self.pdb = """HEADER    First atom from 4NIP (modified)
ATOM      1  N   GLY A   1       4.168   1.038   6.389  1.00 21.86           N
END
"""

    def test_parse_atom_record(self):
        """Test ATOM and HETATM record parsing."""
        reader = pdb_utils.PdbReader()
        f = StringIO(self.pdb)
        f.next()  # skip header
        atom_record = f.next()
        fields = reader.parse_atom_record(atom_record)
        assert fields['record_name'] == 'ATOM'
        assert fields['serial_number'] == 1
        assert fields['atom_name'] == 'N'
        assert fields['alternate_location'] == ''
        assert fields['residue_name'] == 'GLY'
        assert fields['chain'] == 'A'
        assert fields['residue_number'] == 1
        assert fields['insertion_code'] == ''
        assert fields['x'] == 4.168
        assert fields['y'] == 1.038
        assert fields['z'] == 6.389
        assert fields['occupancy'] == 1.
        assert fields['b_factor'] == 21.86
        assert fields['segment'] == ''
        assert fields['element'] == 'N'
        assert fields['charge'] == ''

    def test_pdb_to_pqr(self):
        """Test PDB to PQR conversion."""
        reader = pdb_utils.PdbReader()
        pqr = reader.pdb_to_pqr(StringIO(self.pdb), [-0.35], [2.3])
        ref_pqr = """HEADER    First atom from 4NIP (modified)
ATOM 1 N GLY A 1 4.168 1.038 6.389 -0.35 2.3
END
"""
        assert pqr == ref_pqr, pqr
