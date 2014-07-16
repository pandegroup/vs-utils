"""
AmberTools.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

from collections import OrderedDict
import numpy as np
import os
import shutil
import subprocess
import tempfile

from rdkit import Chem

from pande_gas.utils.pdb_utils import PdbReader


class Antechamber(object):
    """
    Wrapper methods for Antechamber functionality.

    Calculations are carried out in a temporary directory because
    Antechamber writes out several files to disk.

    Parameters
    ----------
    charge_type : str, optional (default 'bcc')
        Antechamber charge type string. Defaults to AM1-BCC charges.
    """
    def __init__(self, charge_type='bcc'):
        self.charge_type = charge_type

        # temporary directory
        self.temp_dir = tempfile.mkdtemp()

    def __del__(self):
        """
        Cleanup.
        """
        shutil.rmtree(self.temp_dir)

    def get_charges_and_radii(self, mol):
        """
        Use Antechamber to calculate partial charges and atomic radii.

        Antechamber requires file inputs and output, so the molecule is
        written to SDF and Antechamber writes out a modified PDB (mpdb)
        containing charge and radius information.

        For multiconformer molecules, charges are calculated using the
        first conformer and then assigned to the remaining conformers.

        Parameters
        ----------
        mol : RDMol
            Molecule.
        """

        # write molecule to temporary file
        _, input_filename = tempfile.mkstemp(suffix='.sdf', dir=self.temp_dir)
        Chem.MolToMolFile(mol, input_filename)

        # calculate charges and radii with Antechamber
        # (antechamber forks and calls SQM also)
        _, output_filename = tempfile.mkstemp(suffix='.mpdb',
                                              dir=self.temp_dir)
        args = ['antechamber', '-i', input_filename, '-fi', 'sdf', '-o',
                output_filename, '-fo', 'mpdb', '-c', self.charge_type]
        subprocess.check_call(args, cwd=self.temp_dir, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

        # extract charges and radii
        reader = ModifiedPdbReader()
        with open(output_filename) as f:
            charges, radii = reader.get_charges_and_radii(f)

        return charges, radii


class PBSA(object):
    """
    Wrapper methods for PBSA functionality.

    Calculations are carried out in a temporary directory because PBSA
    writes out several files to disk.

    Parameters
    ----------
    size : float, optional (default 30.)
        Length of each side of the grid, in Angstroms. Used to calculate
        PBSA parameters xmin, xmax, etc.
    resolution : float, optional (default 0.5)
        Space between grid points, in Angstroms. Corresponds to PBSA space
        parameter.
    nb_cutoff : float, optional (default 5.)
        Cutoff distance for van der Waals interactions. Corresponds to PBSA
        cutnb parameter.
    ionic_strength : float, optional (default 150.)
        Ionic strength of the solvent, in mM. Corresponds to PBSA istrng
        parameter.
    """
    def __init__(self, size=30., resolution=0.5, nb_cutoff=5.,
                 ionic_strength=150.):
        self.size = float(size)
        self.resolution = float(resolution)
        self.nb_cutoff = float(nb_cutoff)
        self.ionic_strength = float(ionic_strength)

        # temporary directory
        self.temp_dir = tempfile.mkdtemp()

    def __del__(self):
        """
        Cleanup.
        """
        shutil.rmtree(self.temp_dir)

    def get_esp_grid(self, pqr):
        """
        Use PBSA to calculate the electrostatic potential grid for a
        molecule.

        PBSA requires file input and output, so we provide a PQR for the
        molecule (one conformer only) and the grid is written is ASCII
        format to pbsa.phi.

        Parameters
        ----------
        pqr : file_like
            Input PQR file.
        """
        # write PQR to disk
        _, pqr_filename = tempfile.mkstemp(suffix='.pqr', dir=self.temp_dir)
        with open(pqr_filename, 'wb') as f:
            f.write(pqr)

        # write PBSA parameter file
        _, param_filename = tempfile.mkstemp(suffix='.in', dir=self.temp_dir)
        with open(param_filename, 'wb') as f:
            f.write(self.get_pbsa_parameter_file())

        # run PBSA
        _, output_filename = tempfile.mkstemp(suffix='.out', dir=self.temp_dir)
        os.remove(output_filename)  # PBSA won't overwrite output filename
        args = ['pbsa', '-i', param_filename, '-o', output_filename, '-pqr',
                pqr_filename]
        try:
            subprocess.check_call(args, cwd=self.temp_dir,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            with open(output_filename) as f:
                output = f.read()
                print output
            raise e

        # extract ESP grid
        with open(os.path.join(self.temp_dir, 'pbsa.phi')) as f:
            grid, center = self.parse_esp_grid(f)

        return grid, center

    def get_pbsa_parameter_file(self):
        """
        Construct a PBSA parameter file.
        """
        params = """
        Calculate ESP for a small molecule
        &cntrl
        inp=0,  ! required for PQR input
        /
        &pb
        npbverb=1,  ! be verbose
        phiout=1, phiform=1,  ! write grid to Amber ASCII file
        istrng={istrng},  ! ionic strength
        space={space},  ! grid spacing
        xmin={xmin}, xmax={xmax},
        ymin={ymin}, ymax={ymax},
        zmin={zmin}, zmax={zmax},
        eneopt=1, cutnb={cutnb},
        /
        """
        delta = self.size / 2.
        params = params.format(
            space=self.resolution,
            istrng=self.ionic_strength,
            xmin=-delta, xmax=delta,
            ymin=-delta, ymax=delta,
            zmin=-delta, zmax=delta,
            cutnb=self.nb_cutoff)

        return params

    def parse_esp_grid(self, grid):
        """
        Parse PBSA ASCII electrostatic potential grid.

        Variables used in the ASCII format:
        * h : grid spacing
        * (gox, goy, goz) : grid origin
        * (xm, ym, zm) : grid dimensions
        * phi : electrostatic potential in kcal/mol-e

        The mapping between one-based grid points (i, j, k) and phi indices
        is p_i = i + xm * (j - 1 + ym * (k - 1)). However, since phi is a
        flattened version of the grid (with Fortran ordering), we can use
        np.reshape to get the 3D grid.

        Spatial coordinates (x, y, z) in the grid are given by
        (gox + h * i, goy + h * j, goz + h * k)

        Parameters
        ----------
        grid : file_like
            Amber ASCII format file.
        """
        h = gox = goy = goz = None
        xm = ym = zm = None
        phi = None
        for line in grid:
            line = line.strip()
            if line.startswith('#'):
                continue
            if h is None:
                h, gox, goy, goz = np.asarray(line.split(), dtype=float)
            elif xm is None:
                xm, ym, zm = np.asarray(line.split(), dtype=int)
            else:
                phi = np.asarray(line.split(), dtype=float)
        grid = np.reshape(phi, (xm, ym, zm), order='F')
        center = (gox, goy, goz)
        print "CENTER", center

        # sanity checks
        assert h == self.resolution

        return grid, center


class ModifiedPdbReader(PdbReader):
    """
    Handle Amber modified PDB files and generate Amber-style PQR files.
    """
    def _parse_atom_record(self, line):
        """
        Parse optional fields in ATOM and HETATM records.

        Amber modified PDB files contain charge, radius and atom type
        information in the fields following the x, y, z coordinates for
        atoms.

        Parameters
        ----------
        line : str
            Amber modified PDB ATOM or HETATM line.
        """
        fields = OrderedDict()
        charge, radius, amber_type = line[54:].strip().split()
        fields['charge'] = charge
        fields['radius'] = radius
        fields['amber_type'] = amber_type

        return fields

    def get_charges_and_radii(self, mpdb):
        """
        Extract atomic charges and radii from an Antechamber modified PDB
        file.

        Parameters
        ----------
        mpdb : file_like
            Antechamber modified PDB file.
        """
        charges = []
        radii = []
        for line in mpdb:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                fields = self.parse_atom_record(line)
                charges.append(fields['charge'])
                radii.append(fields['radius'])
        charges = np.asarray(charges, dtype=float)
        radii = np.asarray(radii, dtype=float)

        return charges, radii
