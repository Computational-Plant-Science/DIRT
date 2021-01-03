from pathlib import Path
from os.path import join

import imageio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import skimage
from scipy import ndimage
import pytest

from thresholding import binary_threshold
from options import DIRTOptions
from segmentation import segment, calculate_histogram, find_marker, find_crown, find_tag, find_tag

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
def test_calculate_histogram(test_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', test_file),
        output_directory=join(test_dir, 'output'))
    image = imageio.imread(options.input_file, as_gray=True)
    masked = binary_threshold(image.astype(np.uint8))
    labeled, _ = ndimage.label(masked) # extract features

    # calculate histogram
    histogram, bin_edges, comps_x, comps_y = calculate_histogram(labeled)

    # for visual sanity check
    plt.bar(bin_edges[:-1], histogram, width=np.shape(labeled)[1] / 100)
    plt.savefig(join(test_dir, 'output', f"{options.input_stem}.histogram.png"))
    plt.clf()


@pytest.mark.parametrize("test_file", test_files)
def test_find_marker(test_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', test_file),
        output_directory=join(test_dir, 'output'))
    image = imageio.imread(options.input_file, as_gray=True)
    masked = binary_threshold(image.astype(np.uint8))
    labeled, _ = ndimage.label(masked)
    histogram, bin_edges, comps_x, comps_y = calculate_histogram(labeled)

    # find marker
    marker, marker_ratio, marker_left, marker_right, marker_bottom, marker_top, _ = find_marker(
        labeled.copy(),
        histogram,
        comps_x,
        comps_y)

    # for visual sanity check, overlay bounding box on original
    fig, ax = plt.subplots(1)
    ax.imshow(image)
    rect = patches.Rectangle((marker_top, marker_left), marker_bottom - marker_top, marker_right - marker_left, edgecolor='r')
    rect.set_alpha(0.4)
    ax.add_patch(rect)
    plt.savefig(join(test_dir, 'output', f"{options.input_stem}.marker.bounding.png"))
    plt.clf()


@pytest.mark.parametrize("test_file", test_files)
def test_find_tag(test_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', test_file),
        output_directory=join(test_dir, 'output'))
    image = imageio.imread(options.input_file, as_gray=True)
    masked = binary_threshold(image.astype(np.uint8))
    labeled, _ = ndimage.label(masked)
    histogram, bin_edges, comps_x, comps_y = calculate_histogram(labeled)
    crown, crown_list, crown_left, crown_right, crown_bottom, crown_top = find_crown(
        labeled.copy(),
        histogram)
    circle, circle_ratio, circle_max, circle_min, circle_top, circle_bottom, _ = find_marker(
        labeled.copy(),
        histogram,
        comps_x,
        comps_y)
    circle_width = float(circle_max) - float(circle_min)
    circle_height = float(circle_bottom) - float(circle_top)
    tag, tag_list, tag_left, tag_right, tag_bottom, tag_top = find_tag(
        labeled.copy(),
        histogram,
        [circle, crown],
        crown,
        crown_list)

    # for visual sanity check, overlay bounding box on original
    fig, ax = plt.subplots(1)
    ax.imshow(image)
    tag_patch = patches.Rectangle((tag_top, tag_left), tag_bottom - tag_top, tag_right - tag_left, edgecolor='r')
    tag_patch.set_alpha(0.4)
    ax.add_patch(tag_patch)
    plt.savefig(join(test_dir, 'output', f"{options.input_stem}.tag.bounding.png"))
    plt.clf()


@pytest.mark.parametrize("test_file", test_files)
def test_find_crown(test_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', test_file),
        output_directory=join(test_dir, 'output'))
    image = imageio.imread(options.input_file, as_gray=True)
    masked = binary_threshold(image.astype(np.uint8))
    labeled, _ = ndimage.label(masked)
    histogram, bin_edges, comps_x, comps_y = calculate_histogram(labeled)

    # find crown
    crown, crown_list, crown_left, crown_right, crown_bottom, crown_top = find_crown(
        labeled.copy(),
        histogram)

    # for visual sanity check, overlay bounding box on original
    fig, ax = plt.subplots(1)
    ax.imshow(image)
    rect = patches.Rectangle((crown_top, crown_left), crown_bottom - crown_top, crown_right - crown_left, edgecolor='r')
    rect.set_alpha(0.4)
    ax.add_patch(rect)
    plt.savefig(join(test_dir, 'output', f"{options.input_stem}.crown.bounding.png"))
    plt.clf()


@pytest.mark.parametrize("test_file", test_files)
@pytest.mark.skip("until fixed")
def test_correct_for_stem(test_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', test_file),
        output_directory=join(test_dir, 'output'))
    image = imageio.imread(options.input_file, as_gray=True)
    masked = binary_threshold(image.astype(np.uint8))
    labeled, _ = ndimage.label(masked)
    histogram, bin_edges, comps_x, comps_y = calculate_histogram(labeled)
    root_index, root_list, crown_left, crown_right, crown_bottom, crown_top = find_crown(labeled.copy(), histogram)

    # correct for stem
    corrected, histogram, stem_left, stem_right, stem_bottom, stem_top = find_tag(
        labeled.copy(),
        histogram,
        [-1, -1, root_index],
        crown_left,
        crown_right,
        crown_bottom,
        crown_top,
        root_index,
        root_list)

    # for visual sanity check, overlay bounding box on original
    fig, ax = plt.subplots(1)
    ax.imshow(image)
    rect = patches.Rectangle((stem_top, stem_left), stem_bottom - stem_top, stem_right - stem_left, edgecolor='r')
    rect.set_alpha(0.4)
    ax.add_patch(rect)
    plt.savefig(join(test_dir, 'output', f"{options.input_stem}.stem.corrected.png"))
    plt.clf()

    # for visual sanity check
    # imageio.imwrite(join(test_dir, 'output', f"{options.input_stem}.stem.corrected.png"), skimage.img_as_uint(corrected))


@pytest.mark.parametrize("test_file", test_files)
def test_preprocess_image(test_file):
    options = DIRTOptions(
        input_file=join(test_dir, 'data', test_file),
        output_directory=join(test_dir, 'output'))
    results = segment(options)
    print(results)
