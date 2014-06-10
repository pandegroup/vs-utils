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
    pixels = np.asarray(image, dtype='uint8')
    #pixels = np.atleast_3d(pixels)  # maintain channels convention
    return pixels


def downscale(image, max_size):
    """
    Shrink an image while maintaining aspect ratio.

    Parameters
    ----------
    image : Image
        Image to rescale.
    max_size : int
        Maximum image size in any dimension.
    """
    image.thumbnail((max_size, max_size), resample=Image.ANTIALIAS)


def pad(pixels, shape, fill=255):
    """
    Pad an image, where the first two axes are height and width,
    respectively.

    Parameters
    ----------
    pixels : ndarray
        Array of pixels, possibly with multiple channels (on third axis).
    shape : tuple of ints
        Desired height and width.
    fill : int
        Intensity value for added pixels.
    """
    current = pixels.shape[:2]
    pad_width = []
    for i in xrange(2):
        diff = shape[i] - current[i]
        n, m = divmod(diff, 2)
        m += n
        pad_width.append((n, m))
    while len(pad_width) < pixels.ndim:
        pad_width.append((0, 0))
    padded = np.pad(pixels, pad_width, 'constant', constant_values=fill)
    return padded


def pad_batch(images, fill=255):
    """
    Pad a set of images to the minimum common shape.

    Parameters
    ----------
    images : list of ndarrays
        Images to pad.
    fill : int
        Intensity value for added pixels.
    """
    h = 0
    w = 0
    for image in images:
        h = max(h, image.shape[0])
        w = max(w, image.shape[1])
    shape = (len(images), h, w)
    if images[0].ndim > 2:
        shape += tuple(images[0].shape[2:])
    padded = np.zeros(shape, dtype='uint8')
    for i, image in enumerate(images):
        padded[i] = pad(image, (h, w), fill)
    if padded.ndim == 3:
        padded = np.expand_dims(padded, 3)
    return padded
