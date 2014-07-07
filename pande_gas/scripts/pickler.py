"""
Convenience script to convert featurized hdf5 files into pickled files of numpy arrays for easier processing.
"""

import argparse
import cPickle as pickle
from featurize import HelpFormatter
import h5py

def main():
  parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
  parser.add_argument('input', help='Input featurized (HDF5) file.')
  parser.add_argument('output', help='Output filename (no extension).')

  args = parser.parse_args()
  f = h5py.File(args.input, "r")

  # Extract compound names
  names = f['names'][:]
  features = f['features'][:]

  with open(args.output + "_names.p", "wb") as name_file:
    pickle.dump(names, name_file)
  with open(args.output + "_features.p", "wb") as feature_file:
    pickle.dump(features, feature_file)

if __name__ == '__main__':
  main()
