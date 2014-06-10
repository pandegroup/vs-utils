"""
Test image featurizer.
"""
from rdkit import Chem

from pande_gas.features import images


def test_images():
    """Test MolImage."""
    mol = Chem.MolFromSmiles(test_smiles)
    f = images.MolImage(shape=(200, 250))
    rval = f([mol])
    assert rval.shape == (1, 200 * 250 * 3)


def test_images_topo_view():
    """Test MolImage using topo_view."""
    mol = Chem.MolFromSmiles(test_smiles)
    f = images.MolImage(shape=(200, 250), flatten=False)
    rval = f([mol])
    assert rval.shape == (1, 200, 250, 3)

test_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
