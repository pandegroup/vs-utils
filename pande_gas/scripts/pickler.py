"""
Convenience script to convert featurized hdf5 files into pickled files of numpy arrays for easier processing.
"""

import argparse
import cPickle as pickle
from featurize import HelpFormatter
import h5py
import os

def hdf5_to_pickle(inputfile, outputfile):
  f = h5py.File(inputfile, "r")

  # Extract compound names
  names = f['names'][:]
  features = f['features'][:]

  with open(outputfile + "_names.p", "wb") as name_file:
    pickle.dump(names, name_file)
  with open(outputfile + "_features.p", "wb") as feature_file:
    pickle.dump(features, feature_file)

def hdf5_to_pickle_all(directory):
  files = os.listdir(directory)
  h5_files = [f for f in files if f.endswith(".h5")]
  for h5_file in h5_files:
    basename = os.path.splitext(h5_file)[0]
    h5_file = os.path.join(directory, h5_file)
    basename = os.path.join(directory, basename)
    hdf5_to_pickle(h5_file, basename)

def main():
  parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
  parser.add_argument('--input', help='Input featurized (HDF5) file.')
  parser.add_argument('--output', help='Output filename (no extension).')
  parser.add_argument('--directory',
          help='Pickle-ify every hdf5 file in given directory')

  args = parser.parse_args()
  if args.directory != None:
    hdf5_to_pickle_all(args.directory)
  else:
    hdf5_to_pickle(args.input, args.output)

if __name__ == '__main__':
  main()
