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
from segmentation import segment, calculate_histogram, find_marker, find_crown, find_tag, find_tag, find_excised

test_dir = Path(__file__).parent
test_images = [
    'cassava1.jpg',
    'cassava2.jpg',
    'hemp1.jpg',
    'hemp2.jpg',
    'maize1.jpg',
    'maize2.jpg',
    'tomato1.jpg',
    'tomato2.jpg',
]
excised_test_images = [
    'hemp1.jpg',
    'hemp2.jpg',
    'tomato1.jpg',
    'tomato2.jpg',
]


@pytest.mark.parametrize("image_file", test_images)
def test_calculate_histogram(image_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', image_file),
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


@pytest.mark.parametrize("image_file", test_images)
def test_find_marker(image_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', image_file),
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
    plt.savefig(join(test_dir, 'output', f"{options.input_stem}.bounding.marker.png"))
    plt.clf()


@pytest.mark.parametrize("image_file", test_images)
def test_find_tag(image_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', image_file),
        output_directory=join(test_dir, 'output'))
    image = imageio.imread(options.input_file, as_gray=True)
    masked = binary_threshold(image.astype(np.uint8))
    labeled, _ = ndimage.label(masked)
    histogram, bin_edges, comps_x, comps_y = calculate_histogram(labeled)
    crown, crown_list, crown_left, crown_right, crown_bottom, crown_top = find_crown(labeled.copy(), histogram)
    circle, circle_ratio, circle_max, circle_min, circle_top, circle_bottom, _ = find_marker(
        labeled.copy(),
        histogram,
        comps_x,
        comps_y)
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
    plt.savefig(join(test_dir, 'output', f"{options.input_stem}.bounding.tag.png"))
    plt.clf()


@pytest.mark.parametrize("image_file", test_images)
def test_find_crown(image_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', image_file),
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
    plt.savefig(join(test_dir, 'output', f"{options.input_stem}.bounding.crown.png"))
    plt.clf()


@pytest.mark.parametrize("image_file", excised_test_images)
@pytest.mark.skip("until fixed")
def test_find_excised(image_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', image_file),
        output_directory=join(test_dir, 'output'),
        excised_roots=1)
    image = imageio.imread(options.input_file, as_gray=True)
    masked = binary_threshold(image.astype(np.uint8))
    labeled, _ = ndimage.label(masked)
    histogram, bin_edges, comps_x, comps_y = calculate_histogram(labeled)
    crown, crown_list, crown_left, crown_right, crown_bottom, crown_top = find_crown(labeled.copy(), histogram)
    marker, marker_ratio, marker_right, marker_left, marker_top, marker_bottom, _ = find_marker(
        labeled.copy(),
        histogram,
        comps_x,
        comps_y)
    tag, tag_list, tag_left, tag_right, tag_bottom, tag_top = find_tag(
        labeled.copy(),
        histogram,
        [marker, crown],
        crown,
        crown_list)

    # find excised segment
    exRIdx, excised_left, excised_right, excised_bottom, excised_top, imgExRoot, centerPtx, centerPty = find_excised(
        labeled.copy(),
        histogram,
        [crown, marker, tag])

    # for visual sanity check, overlay bounding box on original
    fig, ax = plt.subplots(1)
    ax.imshow(image)
    rect = patches.Rectangle(
        (excised_top, excised_left),
        excised_bottom - excised_top,
        excised_right - excised_left,
        edgecolor='r')
    rect.set_alpha(0.4)
    ax.add_patch(rect)
    plt.savefig(join(test_dir, 'output', f"{options.input_stem}.bounding.excised.png"))
    plt.clf()


@pytest.mark.parametrize("image_file", test_images)
@pytest.mark.skip("until fixed")
def test_correct_for_stem(image_file):
    # setup
    options = DIRTOptions(
        input_file=join(test_dir, 'data', image_file),
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
    plt.savefig(join(test_dir, 'output', f"{options.input_stem}.bounding.stem.png"))
    plt.clf()

    # for visual sanity check
    # imageio.imwrite(join(test_dir, 'output', f"{options.input_stem}.stem.corrected.png"), skimage.img_as_uint(corrected))


@pytest.mark.parametrize("image_file", test_images)
def test_segmentation(image_file):
    options = DIRTOptions(
        input_file=join(test_dir, 'data', image_file),
        output_directory=join(test_dir, 'output'))
    image = imageio.imread(options.input_file, as_gray=True)
    masked = binary_threshold(image.astype(np.uint8))
    results = segment(options, masked)
    print(results)
