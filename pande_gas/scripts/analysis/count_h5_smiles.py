import h5py
import sys

"""
Count unique and total SMILES
"""

smiles = set()
total = 0
for filename in sys.argv[1:]:
  with h5py.File(filename, 'r') as f:
    smiles.update(f['smiles'][:])
    total += f['smiles'].shape[0]
print len(smiles), total
