"""
Test topological fingerprints.
"""
from rdkit import Chem

from pande_gas.features import fingerprints as fp


def test_circular_fingerprints():
    """Test CircularFingerprint."""
    mol = Chem.MolFromSmiles(test_smiles)
    f = fp.CircularFingerprint()
    rval = f([mol])
    assert rval.shape == (1, f.size)

test_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
