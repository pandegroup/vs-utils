#!/usr/bin/env python

import argparse

from pande_gas.utils.molecule_net import AssayXMLParser


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
                        help='Input XML file(s).')
    return parser.parse_args(input_args)


def main(filenames):
    """
    Parse XML assay descriptions.

    Parameters
    ----------
    filenames : list
        Filenames.
    """
    for filename in filenames:
        parser = AssayXMLParser(filename)

if __name__ == '__main__':
    args = parse_args()
    main(args.input)
