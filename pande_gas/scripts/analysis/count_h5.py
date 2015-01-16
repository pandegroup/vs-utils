import h5py
import numpy as np
import sys

"""
Count actives and decoys
"""
count = []
actives = []

filenames = sys.argv[2:]  # prefix is sys.argv[1]
for filename in sys.argv[2:]:
  name = filename.split('-circular.h5')[0]
  with h5py.File(filename, 'r') as f:
    try:
      int(name)
      name = 'aid' + name
    except ValueError:
      pass
    assert f['X'].shape[0] == f['smiles'].shape[0] == f['y'].shape[0]
    a = np.count_nonzero(f['y'])
    d = f['y'].shape[0] - a
    count.append(f['y'].shape[0])
    actives.append(a)
    name = sys.argv[1] + '-' + name
    print '{}\t{}\t{}'.format(name, a, d)

print 'DATASETS', len(count)
count = np.asarray(count, dtype=int)
actives = np.asarray(actives, dtype=int)
print 'COUNT', np.mean(count), np.std(count)
rate = np.true_divide(actives, count)
print 'ACT', np.mean(rate), np.std(rate)
