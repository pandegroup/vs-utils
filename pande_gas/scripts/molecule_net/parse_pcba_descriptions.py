#!/usr/bin/env python
"""
Parse PCBA assay descriptions.
"""
import argparse
import glob

from pande_gas.utils.molecule_net import PcbaJsonParser, PcbaXmlParser


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
  return parser.parse_args(input_args)


def main(filenames, input_format='json'):
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
      parser = PcbaJsonParser(filename)
      data = parser.root
      # TODO(rbharath): Swap these out for parser method calls.
      try:
        print
        print "###########################################################"
        print
        if "comment" in data:
          print "comment: " + str(data["comment"])
          print
        if "xref" in data:
          print "xref: " + str(data["comment"])
          print
        if "name" in data:
          print "name: " + str(parser.get_name())
          print
        if "aid_source" in data:
          print "aid_source: " + str(data["aid_source"])
          print
        if "results" in data:
          print "results: " + str(data["results"])
          print
          for entry in data["results"]:
            print "  results -- " + entry["name"]
        if "aid" in data:
          print "aid: " + str(parser.get_aid())
          print
        if "revision" in data:
          print "revision: " + str(data["revision"])
          print
        if "activity_outcome_method" in data:
          print "activity_outcome_method: " + str(data["activity_outcome_method"])
          print
        if "description" in data:
          print "description: " + str(parser.get_description())
          print
        print "###########################################################"
        print
      except:
        print "Exception in parsing..."
        print "###########################################################"
    elif input_format == 'xml':
      parser = PcbaXmlParser(filename)
    else:
      raise NotImplementedError(
          'Unrecognized input format "{}"'.format(input_format))

if __name__ == '__main__':
  args = parse_args()
  print args.input
  if type(args.input) is str:
    filenames = glob.glob(args.input)
  elif type(args.input) is list:
    filenames = args.input
  else:
    raise ValueError("--input must be list or string!")
  main(filenames, args.format)
