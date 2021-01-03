from os.path import join
from pathlib import Path

import imageio
import numpy as np
import skimage
import cv2
import pytest

from thresholding import binary_threshold, simple_threshold, adaptive_threshold_mean, adaptive_threshold_gaussian, otsu_threshold
from options import DIRTOptions

test_dir = Path(__file__).parent
test_files = [
    'cassava1.jpg',
    'cassava2.jpg',
    'hemp1.jpg',
    'hemp2.jpg',
    'maize1.jpg',
    'maize2.jpg',
    'tomato1.jpg',
    'tomato2.jpg',
]


@pytest.mark.parametrize("test_file", test_files)
@pytest.mark.skip(reason="unused")
def test_simple_threshold(test_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', test_file),
        output_directory=join(test_dir, 'output'))
    image = imageio.imread(options.input_file, as_gray=True)

    # apply threshold
    thresholded = simple_threshold(image, 90)

    # make sure image sizes are equal
    assert np.shape(image)[0] == np.shape(thresholded)[0]
    assert np.shape(image)[1] == np.shape(thresholded)[1]

    # for visual sanity check
    imageio.imwrite(join(test_dir, 'output', f"{options.input_stem}.simple.threshold.png"), thresholded)


@pytest.mark.parametrize("test_file", test_files)
@pytest.mark.skip(reason="unused")
def test_mean_adaptive_threshold(test_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', test_file),
        output_directory=join(test_dir, 'output'))
    image = cv2.imread(options.input_file, cv2.IMREAD_GRAYSCALE)

    # apply threshold
    thresholded = adaptive_threshold_mean(image)

    # make sure image sizes are equal
    assert np.shape(image)[0] == np.shape(thresholded)[0]
    assert np.shape(image)[1] == np.shape(thresholded)[1]

    # for visual sanity check
    imageio.imwrite(join(test_dir, 'output', f"{options.input_stem}.adaptive.threshold.mean.png"), thresholded)


@pytest.mark.parametrize("test_file", test_files)
@pytest.mark.skip(reason="unused")
def test_gaussian_adaptive_threshold(test_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', test_file),
        output_directory=join(test_dir, 'output'))
    image = cv2.imread(options.input_file, cv2.IMREAD_GRAYSCALE)

    # apply threshold
    thresholded = adaptive_threshold_gaussian(image)

    # make sure image sizes are equal
    assert np.shape(image)[0] == np.shape(thresholded)[0]
    assert np.shape(image)[1] == np.shape(thresholded)[1]

    # for visual sanity check
    imageio.imwrite(join(test_dir, 'output', f"{options.input_stem}.adaptive.threshold.gaussian.png"), thresholded)


@pytest.mark.parametrize("test_file", test_files)
def test_otsu_threshold(test_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', test_file),
        output_directory=join(test_dir, 'output'))
    image = cv2.imread(options.input_file, cv2.IMREAD_GRAYSCALE)

    # apply threshold
    thresholded = otsu_threshold(image)

    # make sure image sizes are equal
    assert np.shape(image)[0] == np.shape(thresholded)[0]
    assert np.shape(image)[1] == np.shape(thresholded)[1]

    # for visual sanity check
    imageio.imwrite(join(test_dir, 'output', f"{options.input_stem}.otsu.threshold.png"), thresholded)


@pytest.mark.parametrize("test_file", test_files)
def test_binary_threshold(test_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', test_file),
        output_directory=join(test_dir, 'output'))
    image = imageio.imread(options.input_file, as_gray=True)

    # apply threshold
    masked = binary_threshold(image.astype(np.uint8))

    # for visual sanity check
    imageio.imwrite(join(test_dir, 'output', f"{options.input_stem}.masked.png"), skimage.img_as_uint(masked))
