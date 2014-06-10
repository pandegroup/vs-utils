"""
Molecular depiction features.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import io
import numpy as np
from PIL import Image


def load(string):
    """
    Load an image from a file or binary string.

    Parameters
    ----------
    string : str
        Filename or binary string.
    """
    try:
        im = Image.open(string)
    except TypeError:
        b = io.BytesIO(string)
        im = Image.open(b)
    return im


def get_pixels(image, mode=None):
    """
    Extract pixels from an image, possibly after converting to the given
    mode.

    Parameters
    ----------
    image : PIL Image
        Image.
    mode : str, optional
        Image mode. For example, 'RGB' or 'P' (8-bit).
    """
    if mode is not None and image.mode != mode:
        image = image.convert(mode)
    pixels = np.asarray(image)
    return pixels


def downscale(image, max_size):
    """
    Shrink an image while maintaining aspect ratio. Returns a copy of the
    original image.

    Parameters
    ----------
    image : Image
        Image to rescale.
    max_size : int
        Maximum image size in any dimension.
    """
    if max(image.size) <= max_size:
        return image.copy()
    im = image.copy()
    im.thumbnail((max_size, max_size), resample=Image.ANTIALIAS)
    return im


def pad(image, shape, fill=255):
    """
    Pad an image, where the first two axes are height and width,
    respectively. Returns a copy of the original image.

    Parameters
    ----------
    image : PIL Image
        Image.
    shape : tuple of ints
        Desired height and width.
    fill : int, optional (default 255)
        Intensity value for added pixels.
    """
    pixels = get_pixels(image)
    current = pixels.shape[:2]
    assert current[0] <= shape[0] and current[1] <= shape[1]
    pad_width = []
    for i in xrange(2):
        diff = shape[i] - current[i]
        n, m = divmod(diff, 2)
        m += n
        pad_width.append((n, m))
    while len(pad_width) < pixels.ndim:
        pad_width.append((0, 0))
    padded = np.pad(pixels, pad_width, 'constant', constant_values=fill)
    im = Image.fromarray(padded, mode=image.mode)
    return im
