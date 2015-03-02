"""
Test featurizer class.
"""
import numpy as np
import unittest

from rdkit import Chem

from vs_utils.features import MolPreparator
from vs_utils.features.basic import MolecularWeight
from vs_utils.utils.parallel_utils import LocalCluster
from vs_utils.utils.rdkit_utils import conformers


class TestFeaturizer(unittest.TestCase):
    """
    Tests for Featurizer.
    """
    def setUp(self):
        """
        Set up tests.
        """
        smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
        mol = Chem.MolFromSmiles(smiles)
        engine = conformers.ConformerGenerator(max_conformers=1)
        self.mol = engine.generate_conformers(mol)
        assert self.mol.GetNumConformers() > 0

    def test_featurizer(self):
        """
        Test basic functionality of Featurizer.
        """
        f = MolecularWeight()
        rval = f([self.mol])
        assert rval.shape == (1, 1)

    def test_flatten_conformers(self):
        """
        Calculate molecule-level features for a multiconformer molecule.
        """
        f = MolecularWeight()
        rval = f([self.mol])
        assert rval.shape == (1, 1)

    def test_parallel(self):
        """
        Test parallel featurization.
        """
        cluster = LocalCluster(1)
        f = MolecularWeight()
        rval = f([self.mol])
        parallel_rval = f([self.mol], parallel=True,
                          client_kwargs={'cluster_id': cluster.cluster_id})
        assert np.array_equal(rval, parallel_rval)


class TestMolPreparator(unittest.TestCase):
    """
    Test MolPreparator.
    """
    def setUp(self):
        """
        Set up tests.
        """
        smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
        mol = Chem.MolFromSmiles(smiles)
        engine = conformers.ConformerGenerator(max_conformers=1)
        self.mol = engine.generate_conformers(mol)
        assert self.mol.GetNumConformers() > 0
        self.preparator = MolPreparator()

    def test_identity(self):
        """
        Test MolPreparator does nothing by default.
        """
        mol = self.preparator(self.mol)
        assert mol.GetNumAtoms() == self.mol.GetNumAtoms()
        assert mol.GetNumConformers() == self.mol.GetNumConformers()

    def test_ionize(self):
        """
        Test MolPreparator with ionize=True.
        """
        self.preparator.set_ionize(True)
        mol = self.preparator(self.mol)
        mol_charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
        ref_charge = sum([atom.GetFormalCharge()
                          for atom in self.mol.GetAtoms()])
        assert mol_charge != ref_charge

    def test_align(self):
        """
        Test MolPreparator with align=True.

        Really just make sure that coordinates are sane, since the alignment
        tends to flatten molecules if there are hydrogens present.
        """
        def _check_flat(this_mol):
            """
            Check whether a molecule is flat.

            Parameters
            ----------
            this_mol : RDKit Mol
                Molecule to check.
            """
            for conf in this_mol.GetConformers():
                coords = []
                for atom in this_mol.GetAtoms():
                    coords.append(list(conf.GetAtomPosition(atom.GetIdx())))
                coords = np.asarray(coords, dtype=float)
                for i in xrange(coords.shape[1]):
                    assert not np.allclose(coords[:, i], 0)

        # test without hydrogens
        self.preparator.set_align(True)
        mol = self.preparator(self.mol)
        _check_flat(mol)

        # test with hydrogens
        self.preparator.set_add_hydrogens(True)
        mol = self.preparator(self.mol)
        _check_flat(mol)

    def test_add_hydrogens(self):
        """
        Test MolPreparator with add_hydrogens=True.
        """
        self.preparator.set_add_hydrogens(True)
        ref_mol = Chem.RemoveHs(self.mol)
        mol = self.preparator(ref_mol)
        assert mol.GetNumAtoms() > ref_mol.GetNumAtoms()
