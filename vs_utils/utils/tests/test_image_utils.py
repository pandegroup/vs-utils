"""
Tests for image utilities.
"""
import numpy as np
from PIL import Image
import tempfile
import unittest

from vs_utils.utils import image_utils


class TestImageUtils(unittest.TestCase):
    """
    Test image_utils.
    """
    def setUp(self):
        """
        Set up tests.
        """
        pixels = np.random.randint(255, size=(5, 5, 3))
        self.pixels = np.asarray(pixels, dtype='uint8')

    def test_get_pixels(self):
        """
        Read pixels from image.
        """
        im = Image.fromarray(self.pixels, mode='RGB')
        assert np.array_equal(self.pixels, image_utils.get_pixels(im))

    def test_load_file(self):
        """
        Load an image from a file.
        """
        im = Image.fromarray(self.pixels, mode='RGB')
        _, filename = tempfile.mkstemp(suffix='.png')
        im.save(filename)
        im = image_utils.load(filename)
        assert np.array_equal(self.pixels, image_utils.get_pixels(im))

    def test_load_string(self):
        """
        Load an image from binary string.
        """
        im = Image.fromarray(self.pixels, mode='RGB')
        _, filename = tempfile.mkstemp(suffix='.png')
        im.save(filename)
        with open(filename) as f:
            string = f.read()
        im = image_utils.load(string)
        assert np.array_equal(self.pixels, image_utils.get_pixels(im))

    def test_downscale(self):
        """
        Downscale image while maintaining aspect ratio.
        """
        pixels = np.random.randint(255, size=(20, 16, 3))
        pixels = np.asarray(pixels, dtype='uint8')
        im = Image.fromarray(pixels, mode='RGB')
        im = image_utils.downscale(im, 10)
        assert im.size == (8, 10)

    def test_pad(self):
        """Pad an image."""
        im = Image.fromarray(self.pixels, mode='RGB')
        im = image_utils.pad(im, (7, 8))
        assert im.size == (8, 7)

    def test_pad_fail(self):
        """
        Attempt to pad an image with desired size smaller than original size.
        """
        im = Image.fromarray(self.pixels, mode='RGB')
        try:
            image_utils.pad(im, (4, 3))
        except AssertionError:
            return True
        raise AssertionError
