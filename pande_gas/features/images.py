"""
Molecule images.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import io
from PIL import Image
import os
import subprocess


def image_from_smiles(smiles):
    """
    Generate a PNG image from a SMILES string. Uses OpenBabel to generate
    an SVG (which gives uniform scaling) and ImageMagick convert to
    generate a PNG.

    Note that we call OpenBabel using Popen to avoid entanglements with the
    GPL.

    See http://stackoverflow.com/questions/13332268 for details on the
    echo function.

    Parameters
    ----------
    smiles : str
        Canonical SMILES string.
    """
    devnull = open(os.devnull, 'w')
    echo = lambda string: subprocess.Popen(['echo', string],
                                           stdout=subprocess.PIPE).stdout
    svg_args = ['obabel', '-ican', '-osvg', '-d', '-xd']
    svg = subprocess.check_output(svg_args, stdin=echo(smiles), stderr=devnull)
    png_args = ['convert', '-alpha', 'off', 'svg:-', 'png:-']
    png = subprocess.check_output(png_args, stdin=echo(svg), stderr=devnull)
    b = io.BytesIO(png)
    im = Image.open(b)
    return im
