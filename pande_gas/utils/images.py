"""
Molecular depiction features.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import io
import numpy as np
from PIL import Image


def read_pixels(filename, mode=None, max_size=None):
    """
    Read pixels from an image file.

    Parameters
    ----------
    filename : str
        Image filename.
    """
    im = Image.open(filename)
    return get_pixels(im, mode, max_size)


def read_pixels_from_string(string, mode=None, max_size=None):
    """Read pixels from a binary string.

    Parameters
    ----------
    string : str
        Binary image string.
    """
    b = io.BytesIO(string)
    im = Image.open(b)
    return get_pixels(im, mode, max_size)


def get_pixels(image, mode=None, max_size=None):
    """
    Extract pixels from an image.

    Parameters
    ----------
    image : PIL Image
        Image.
    mode : str or None
        Image mode. For example, 'RGB' or 'P' (8-bit).
    max_size : int or None
        Scale images such that no dimension is larger than max_size.
    """
    if mode is not None:
        if image.mode != mode:
            image = image.convert(mode)
    if max_size is not None:
        downscale(image, max_size)
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


def pad_batch(images, fill=255):
    """
    Pad a set of images to the minimum common shape.

    Parameters
    ----------
    images : list
        Images to pad.
    fill : int, optional (default 255)
        Intensity value for added pixels.
    """
    h = 0
    w = 0
    for image in images:
        h = max(h, image.size[1])
        w = max(w, image.size[0])
    padded = []
    for i, image in enumerate(images):
        padded.append(pad(image, (h, w), fill))
    return padded
