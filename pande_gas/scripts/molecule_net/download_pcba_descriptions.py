#!/usr/bin/env python
"""
Download PCBA assay descriptions.
"""
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
    parser.add_argument('-f', '--format', choices=['json', 'xml'])
    return parser.parse_args(input_args)


def main(filename, output_format='json'):
    """
    Download PCBA JSON descriptions.

    Parameters
    ----------
    filename : str
        Filename containing AIDs.
    output_format : str, optional (default 'json')
        Output file format.
    """
    with open(filename) as f:
        aids = [int(line) for line in f]
    engine = PubChem()
    for aid in aids:
        data = engine.get_assay_description(aid, output_format)
        with open('aid{}.{}'.format(aid, output_format), 'wb') as f:
            f.write(data)

if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.format)
