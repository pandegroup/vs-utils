"""
Handle PDB files.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

from collections import OrderedDict


class PdbReader(object):
    """
    Handle PDB files.

    Also supports conversion from PDB to Amber-style PQR files.
    """
    def parse_atom_record(self, line):
        """
        Extract fields from a PDB ATOM or HETATM record.

        See http://deposit.rcsb.org/adit/docs/pdb_atom_format.html.

        Parameters
        ----------
        line : str
            PDB ATOM or HETATM line.
        """
        assert line.startswith('ATOM') or line.startswith('HETATM')

        fields = OrderedDict()
        fields['record_name'] = line[:6]
        fields['serial_number'] = int(line[6:11])
        fields['atom_name'] = line[12:16]
        fields['alternate_location'] = line[16]
        fields['residue_name'] = line[17:20]
        fields['chain'] = line[21]
        fields['residue_number'] = int(line[22:26])
        fields['insertion_code'] = line[26]
        fields['x'] = float(line[30:38])
        fields['y'] = float(line[38:46])
        fields['z'] = float(line[46:54])

        # parse additional fields
        fields.update(self._parse_atom_record(line))

        # strip extra whitespace from fields
        for key in fields.keys():
            try:
                fields[key] = fields[key].strip()
            except AttributeError:
                pass

        return fields

    def _parse_atom_record(self, line):
        """
        Parse optional fields in ATOM and HETATM records.

        Parameters
        ----------
        line : str
            PDB ATOM or HETATM line.
        """
        fields = OrderedDict()
        try:
            fields['occupancy'] = float(line[54:60])
            fields['b_factor'] = float(line[60:66])
            fields['segment'] = line[72:76]
            fields['element'] = line[76:78]
            fields['charge'] = line[78:80]
        except IndexError:
            pass

        return fields

    def pdb_to_pqr(self, pdb, charges, radii):
        """
        Convert PDB to Amber-style PQR by adding charge and radius
        information. See p. 68 of the Amber 14 Reference Manual.

        Parameters
        ----------
        pdb : file_like
            PDB file.
        charges : array_like
            Atomic partial charges.
        radii : array_like
            Atomic radii.
        """

        # only certain PDB fields are used in the Amber PQR format
        pdb_fields = ['record_name', 'serial_number', 'atom_name',
                      'residue_name', 'chain', 'residue_number', 'x', 'y', 'z']

        i = 0
        pqr = ''
        for line in pdb:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                fields = self.parse_atom_record(line)

                # charge and radius are added after x, y, z coordinates
                pqr_fields = []
                for field in pdb_fields:
                    value = fields[field]
                    if value == '':
                        value = '?'
                    pqr_fields.append(str(value))
                pqr_fields.append(str(charges[i]))
                pqr_fields.append(str(radii[i]))
                line = ' '.join(pqr_fields) + '\n'
                i += 1  # update atom count
            pqr += line

        # check that we covered all the atoms
        assert i == len(charges) == len(radii)

        return pqr
