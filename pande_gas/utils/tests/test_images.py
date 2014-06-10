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


def test_read_pixels():
    """Read pixels from image file."""
    p = np.random.randint(255, size=(5, 5, 3))
    p = np.asarray(p, dtype='uint8')
    im = Image.fromarray(p, mode='RGB')
    _, filename = tempfile.mkstemp(suffix='.png')
    im.save(filename)
    assert np.array_equal(p, images.read_pixels(filename))


def test_read_pixels_from_string():
    """Read pixels from binary string."""
    p = np.random.randint(255, size=(5, 5, 3))
    p = np.asarray(p, dtype='uint8')
    im = Image.fromarray(p, mode='RGB')
    _, filename = tempfile.mkstemp(suffix='.png')
    im.save(filename)
    with open(filename) as f:
        string = f.read()
    assert np.array_equal(p, images.read_pixels_from_string(string))


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


def test_pad_batch():
    """Pad a set of images to a common shape."""
    imgs = []
    h = 0
    w = 0
    for i in xrange(10):
        size = tuple(np.random.randint(5, 10, size=2)) + (3,)
        if size[0] > h:
            h = size[0]
        if size[1] > w:
            w = size[1]
        p = np.random.randint(255, size=size)
        p = np.asarray(p, dtype='uint8')
        im = Image.fromarray(p, mode='RGB')
        imgs.append(im)
    imgs = images.pad_batch(imgs)
    for im in imgs:
        assert im.size == (w, h)
