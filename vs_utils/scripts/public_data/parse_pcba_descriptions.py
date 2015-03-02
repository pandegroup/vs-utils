#!/usr/bin/env python
"""
Parse PCBA assay descriptions.
"""
import argparse
import glob

from pande_gas.utils.public_data import PcbaJsonParser, PcbaPandasHandler


def parse_args(input_args=None):
  """
  Parse command-line arguments.

  Parameters
  ----------
  input_args : list, optional
      Input arguments. If not provided, defaults to sys.argv[1:].
  """
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input', required=1, nargs='+',
                      help='Input file(s) containing assay description(s). '
                           'Can be in glob format.')
  parser.add_argument('-f', '--format', choices=['json', 'xml'],
                      default='json',
                      help='Input file format.')
  parser.add_argument("--out", help="Location of CSV output.",
                      required=True)
  return parser.parse_args(input_args)


def main(filenames, input_format='json', outfile="~/out.txt"):
  """
  Parse PCBA assay descriptions.

  Parameters
  ----------
  filenames : list
      Filenames containing assay descriptions.
  input_format : str, optional (default 'json')
      Input file format.
  """
  for filename in filenames:
    if input_format == 'json':
      handler = PcbaPandasHandler()
      handler.add_dataset(filename)
    else:
      raise NotImplementedError(
          'Unrecognized input format "{}"'.format(input_format))
  handler.write(out_file)

if __name__ == '__main__':
  args = parse_args()
  print args.input
  if type(args.input) is str:
    filenames = glob.glob(args.input)
  elif type(args.input) is list:
    filenames = args.input
  else:
    raise ValueError("--input must be list or string!")
  main(filenames, args.format, args.out)
