"""
Test feature-based datasets.
"""
import tempfile

from pylearn2.config import yaml_parse

from pande_gas.utils import rdkit_utils as rd


def test_featurizer_dataset():
    """Test FeaturizerDataset with Train."""
    _, filename = tempfile.mkstemp(suffix='.smi')
    with open(filename, 'wb') as f:
        f.write(test_smiles)
    trainer = yaml_parse.load(test_featurizer_dataset_yaml %
                              {'filename': filename})
    trainer.main_loop()


def test_featurizer_dataset_conformers():
    """Test FeaturizerDataset with Train and conformer features."""
    _, filename = tempfile.mkstemp(suffix='.smi')
    with open(filename, 'wb') as f:
        f.write(test_smiles)
    mols = rd.read(filename)
    mols = [rd.generate_conformers(mol) for mol in mols]
    _, filename = tempfile.mkstemp(suffix='.sdf')
    rd.write(mols, filename)
    trainer = yaml_parse.load(test_featurizer_dataset_conformers_yaml %
                              {'filename': filename})
    trainer.main_loop()


def test_featurization_dataset_cv():
    """Test FeaturizerDataset with TrainCV."""
    _, filename = tempfile.mkstemp(suffix='.smi')
    with open(filename, 'wb') as f:
        f.write(test_smiles)
    trainer = yaml_parse.load(test_featurizer_dataset_cv_yaml %
                              {'filename': filename})
    trainer.main_loop()


def test_featurizer_dataset_cv_conformers():
    """Test FeaturizerDataset with TrainCV and conformer features."""
    _, filename = tempfile.mkstemp(suffix='.smi')
    with open(filename, 'wb') as f:
        f.write(test_smiles)
    mols = rd.read(filename)
    mols = [rd.generate_conformers(mol) for mol in mols]
    _, filename = tempfile.mkstemp(suffix='.sdf')
    rd.write(mols, filename)
    trainer = yaml_parse.load(test_featurizer_dataset_cv_conformers_yaml %
                              {'filename': filename})
    trainer.main_loop()

test_smiles = r"""CC(=O)OC1=CC=CC=C1C(=O)O aspirin
CN1C(=C(C2=CC=CC=C2S1(=O)=O)O)C(=O)NC3=CC=CC=N3 piroxicam
CC\1=C(C2=C(/C1=C\C3=CC=C(C=C3)S(=O)C)C=CC(=C2)F)CC(=O)O sulindac
C1CN2C(=CC=C2C(=O)C3=CC=CC=C3)C1C(=O)O ketorolac
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O ibuprofen
CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F celecoxib
CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C penicillin_g
"""

test_featurizer_dataset_yaml = """
!obj:pylearn2.train.Train {
    dataset: !obj:pande_gas.features.datasets.FeaturizerDataset {
        mols: %(filename)s,
        featurizers: [
            !obj:pande_gas.features.basic.MolecularWeight {},
        ],
        mol_iterator: 0,
    },
    model: !obj:pylearn2.models.autoencoder.Autoencoder {
        nvis: 1,
        nhid: 2,
        act_enc: 'sigmoid',
        act_dec: 'linear'
    },
    algorithm: !obj:pylearn2.training_algorithms.bgd.BGD {
        batch_size: 2,
        line_search_mode: 'exhaustive',
        conjugate: 1,
        termination_criterion:
            !obj:pylearn2.termination_criteria.EpochCounter {
                    max_epochs: 1,
        },
        cost: !obj:pylearn2.costs.autoencoder.MeanSquaredReconstructionError {
        },
    },
}
"""

test_featurizer_dataset_conformers_yaml = """
!obj:pylearn2.train.Train {
    dataset: !obj:pande_gas.features.datasets.FeaturizerDataset {
        mols: %(filename)s,
        featurizers: [
            !obj:pande_gas.features.basic.MolecularWeight {},
            !obj:pande_gas.features.coulomb_matrices.CoulombMatrix {},
        ],
        mol_iterator: 0,
    },
    model: !obj:pylearn2.models.autoencoder.Autoencoder {
        nvis: 352,
        nhid: 10,
        act_enc: 'sigmoid',
        act_dec: 'linear'
    },
    algorithm: !obj:pylearn2.training_algorithms.bgd.BGD {
        batch_size: 2,
        line_search_mode: 'exhaustive',
        conjugate: 1,
        termination_criterion:
            !obj:pylearn2.termination_criteria.EpochCounter {
                    max_epochs: 1,
        },
        cost: !obj:pylearn2.costs.autoencoder.MeanSquaredReconstructionError {
        },
    },
}
"""

test_featurizer_dataset_cv_yaml = """
!obj:pylearn2.cross_validation.TrainCV {
    dataset_iterator:
        !obj:pylearn2.cross_validation.dataset_iterators.DatasetKFold {
        dataset: !obj:pande_gas.features.datasets.FeaturizerDataset {
            mols: %(filename)s,
            featurizers: [
                !obj:pande_gas.features.basic.MolecularWeight {},
            ],
        },
    },
    model: !obj:pylearn2.models.autoencoder.Autoencoder {
        nvis: 1,
        nhid: 2,
        act_enc: 'sigmoid',
        act_dec: 'linear'
    },
    algorithm: !obj:pylearn2.training_algorithms.bgd.BGD {
        batch_size: 2,
        line_search_mode: 'exhaustive',
        conjugate: 1,
        termination_criterion:
            !obj:pylearn2.termination_criteria.EpochCounter {
                    max_epochs: 1,
        },
        cost: !obj:pylearn2.costs.autoencoder.MeanSquaredReconstructionError {
        },
    },
}
"""

test_featurizer_dataset_cv_conformers_yaml = """
!obj:pylearn2.cross_validation.TrainCV {
    dataset_iterator:
        !obj:pylearn2.cross_validation.dataset_iterators.DatasetKFold {
        dataset: !obj:pande_gas.features.datasets.FeaturizerDataset {
            mols: %(filename)s,
            featurizers: [
                !obj:pande_gas.features.basic.MolecularWeight {},
                !obj:pande_gas.features.coulomb_matrices.CoulombMatrix {},
            ],
        },
    },
    model: !obj:pylearn2.models.autoencoder.Autoencoder {
        nvis: 352,
        nhid: 10,
        act_enc: 'sigmoid',
        act_dec: 'linear'
    },
    algorithm: !obj:pylearn2.training_algorithms.bgd.BGD {
        batch_size: 2,
        line_search_mode: 'exhaustive',
        conjugate: 1,
        termination_criterion:
            !obj:pylearn2.termination_criteria.EpochCounter {
                    max_epochs: 1,
        },
        cost: !obj:pylearn2.costs.autoencoder.MeanSquaredReconstructionError {
        },
    },
}
"""
