"""
Tests for image utilities.
"""
import numpy as np
from PIL import Image
import tempfile

from pande_gas.utils import images


def test_get_pixels():
    """Read pixels from image."""
    p = np.random.randint(255, size=(5, 5, 3))
    p = np.asarray(p, dtype='uint8')
    im = Image.fromarray(p, mode='RGB')
    assert np.array_equal(p, images.get_pixels(im))


def test_load_file():
    """Load an image from a file."""
    p = np.random.randint(255, size=(5, 5, 3))
    p = np.asarray(p, dtype='uint8')
    im = Image.fromarray(p, mode='RGB')
    _, filename = tempfile.mkstemp(suffix='.png')
    im.save(filename)
    im = images.load(filename)
    assert np.array_equal(p, images.get_pixels(im))


def test_load_string():
    """Load an image from binary string."""
    p = np.random.randint(255, size=(5, 5, 3))
    p = np.asarray(p, dtype='uint8')
    im = Image.fromarray(p, mode='RGB')
    _, filename = tempfile.mkstemp(suffix='.png')
    im.save(filename)
    with open(filename) as f:
        string = f.read()
    im = images.load(string)
    assert np.array_equal(p, images.get_pixels(im))


def test_downscale():
    """Downscale image while maintaining aspect ratio."""
    p = np.random.randint(255, size=(20, 16, 3))
    p = np.asarray(p, dtype='uint8')
    im = Image.fromarray(p, mode='RGB')
    im = images.downscale(im, 10)
    assert im.size == (8, 10)


def test_pad():
    """Pad an image."""
    p = np.random.randint(255, size=(5, 5, 3))
    p = np.asarray(p, dtype='uint8')
    im = Image.fromarray(p, mode='RGB')
    im = images.pad(im, (7, 8))
    assert im.size == (8, 7)


def test_pad_fail():
    """
    Attempt to pad an image with desired size smaller than original size.
    """
    p = np.random.randint(255, size=(5, 5, 3))
    p = np.asarray(p, dtype='uint8')
    im = Image.fromarray(p, mode='RGB')
    try:
        images.pad(im, (4, 3))
    except AssertionError:
        return True
    raise AssertionError
