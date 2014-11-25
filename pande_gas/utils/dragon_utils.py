"""
Dragon utilities.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

from cStringIO import StringIO
import numpy as np
import os
import pandas as pd
import subprocess
import tempfile

from pande_gas.utils import SmilesGenerator


class Dragon(object):
    """
    Wrapper for dragon6shell.

    Parameters
    ----------
    subset : str, optional (default '2d')
        Descriptor subset.
    kwargs : dict, optional
        Keyword arguments for SmilesGenerator.
    """
    def __init__(self, subset='2d', **kwargs):
        self.subset = subset
        self.initialized = False
        self.config_filename, self.smiles_engine = None, None
        self.smiles_engine_kwargs = kwargs

    def initialize(self):
        """
        Initialize.

        This is not part of __init__ because it breaks IPython.parallel.
        """
        fd, self.config_filename = tempfile.mkstemp()
        os.close(fd)
        with open(self.config_filename, 'wb') as f:
            f.write(self.get_config())
        self.smiles_engine = SmilesGenerator(**self.smiles_engine_kwargs)
        self.initialized = True

    def __del__(self):
        """
        Cleanup.
        """
        if self.config_filename is not None:
            os.unlink(self.config_filename)

    def get_config(self):
        """
        Get configuration file.
        """
        if self.subset == '2d':
            return """<?xml version="1.0" encoding="utf-8"?>
<DRAGON version="6.0.36" script_version="1" generation_date="2014/11/17">
  <OPTIONS>
    <CheckUpdates value="true"/>
    <SaveLayout value="true"/>
    <ShowWorksheet value="false"/>
    <Decimal_Separator value="."/>
    <Missing_String value="NaN"/>
    <DefaultMolFormat value="1"/>
    <HelpBrowser value="/usr/bin/xdg-open"/>
    <RejectUnusualValence value="false"/>
    <Add2DHydrogens value="false"/>
    <MaxSRforAllCircuit value="19"/>
    <MaxSR value="35"/>
    <MaxSRDetour value="30"/>
    <MaxAtomWalkPath value="2000"/>
    <LogPathWalk value="true"/>
    <LogEdge value="true"/>
    <Weights>
      <weight name="Mass"/>
      <weight name="VdWVolume"/>
      <weight name="Electronegativity"/>
      <weight name="Polarizability"/>
      <weight name="Ionization"/>
      <weight name="I-State"/>
    </Weights>
    <SaveOnlyData value="false"/>
    <SaveLabelsOnSeparateFile value="false"/>
    <SaveFormatBlock value="%b - %n.txt"/>
    <SaveFormatSubBlock value="%b-%s - %n - %m.txt"/>
    <SaveExcludeMisVal value="false"/>
    <SaveExcludeAllMisVal value="false"/>
    <SaveExcludeConst value="false"/>
    <SaveExcludeNearConst value="false"/>
    <SaveExcludeStdDev value="false"/>
    <SaveStdDevThreshold value="0.0001"/>
    <SaveExcludeCorrelated value="false"/>
    <SaveCorrThreshold value="0.95"/>
    <SaveExclusionOptionsToVariables value="false"/>
    <SaveExcludeMisMolecules value="false"/>
    <SaveExcludeRejectedMolecules value="false"/>
  </OPTIONS>
  <DESCRIPTORS>
    <block id="1" SelectAll="true"/>
    <block id="2" SelectAll="true"/>
    <block id="3" SelectAll="true"/>
    <block id="4" SelectAll="true"/>
    <block id="5" SelectAll="true"/>
    <block id="6" SelectAll="true"/>
    <block id="7" SelectAll="true"/>
    <block id="8" SelectAll="true"/>
    <block id="9" SelectAll="true"/>
    <block id="10" SelectAll="true"/>
    <block id="11" SelectAll="true"/>
    <block id="12" SelectAll="true"/>
    <block id="21" SelectAll="true"/>
    <block id="22" SelectAll="true"/>
    <block id="23" SelectAll="true"/>
    <block id="24" SelectAll="true"/>
    <block id="25" SelectAll="true"/>
    <block id="28" SelectAll="true"/>
    <block id="29" SelectAll="true"/>
  </DESCRIPTORS>
  <MOLFILES>
    <molInput value="stdin"/>
    <molInputFormat value="SMILES"/>
  </MOLFILES>
  <OUTPUT>
    <SaveStdOut value="true"/>
    <SaveProject value="false"/>
    <SaveFile value="false"/>
    <logMode value="stderr"/>
  </OUTPUT>
</DRAGON>
"""
        else:
            raise NotImplementedError

    def get_descriptors(self, mols):
        """
        Parameters
        ----------
        mols : array_like
            Molecules.
        """
        if not self.initialized:
            self.initialize()
        smiles = [self.smiles_engine.get_smiles(mol) for mol in mols]
        args = ['dragon6shell', '-s', self.config_filename]
        p = subprocess.Popen(args, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate('\n'.join(smiles))
        if not stdout:
            raise RuntimeError(stderr)
        data, names = self.parse_descriptors(stdout)

        # adjust for skipped molecules
        # descriptors are in same order as smiles
        missing = np.setdiff1d(smiles, names)
        features = []
        i = 0  # index into calculated features
        for this_smiles in smiles:
            if this_smiles in missing:
                features.append(None)
            else:
                assert this_smiles == names[i]  # confirm match
                features.append(data[i])
                i += 1
        assert len(features) == len(mols)
        return np.asarray(features)

    def parse_descriptors(self, string):
        """
        Parse Dragon descriptors.

        Parameters
        ----------
        string : str
            Output from dragon6shell.
        """
        df = pd.read_table(StringIO(string))
        if self.subset == '2d':
            del df['nHBonds'], df['Psi_e_1d'], df['Psi_e_1s']

        # extract names
        names = df['NAME'].values

        # delete No. and NAME columns
        del df['No.'], df['NAME']

        return np.asarray(df, dtype=float), names
