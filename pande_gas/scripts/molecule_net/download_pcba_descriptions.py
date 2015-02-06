#!/usr/bin/env python
"""
Download PCBA assay descriptions.
"""
import os
import argparse

from pubchem_utils import PubChem


def parse_args(input_args=None):
    """
    Parse command-line arguments.

    Parameters
    ----------
    input_args : list, optional
        Input arguments. If not provided, defaults to sys.argv[1:].
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        help='Input file containing AIDs.')
    parser.add_argument('-f', '--format', choices=['json', 'xml'],
                        required=True)
    parser.add_argument('--out',
                        help="Output directory for generated files.",
                        required=True)
    return parser.parse_args(input_args)


def main(filename, output_format, output_dir):
    """
    Download PCBA JSON descriptions.

    Parameters
    ----------
    filename : str
        Filename containing AIDs.
    output_format : str (default 'json')
        Output file format.
    output_dir : str
        Output directory where downloaded files will be written.
    """
    if not os.path.isdir(output_dir):
      raise ValueError("%s is not a valid output directory!" % output_dir)
    with open(filename) as f:
      aids = [int(line) for line in f]
    engine = PubChem()
    for aid in aids:
      out_file = os.path.join(
          output_dir, 'aid{}.{}'.format(aid, output_format))
      data = engine.get_assay_description(aid, output_format)
      with open(out_file, 'wb') as f:
          f.write(data)

if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.format, args.out)
