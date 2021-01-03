import time
from os.path import join
from typing import List

import imageio
import matplotlib.pyplot as plt
import numpy as np
import skimage
from matplotlib import patches
from scipy import ndimage
from scipy.ndimage import distance_transform_edt
from skimage.morphology import medial_axis

from options import DIRTOptions
from results import DIRTResults


def calculate_histogram(image: np.ndarray) -> (np.ndarray, np.array, int, int, List[list], List[list]):
    flattened = image.flatten()
    histogram, bin_edges = np.histogram(flattened, bins=np.max(image) + 1)

    # background can have less pixels than foreground if no markers are in the image
    if len(np.unique(image)) == 2:  # image is binary
        white_pixels = len(np.where(image == 255)[1])
        largest_feature = np.max(histogram)
        if largest_feature == white_pixels:
            largest_feature = np.min(histogram)
        histogram[list(histogram).index(largest_feature)] = 0  # add feature
    else:
        largest_feature = np.max(histogram)
        histogram[list(histogram).index(largest_feature)] = 0

    for i in range(len(histogram)):
        if histogram[i] < 100:
            histogram[i] = 0

    x_components = [[] for _ in range(len(histogram))]
    y_components = [[] for _ in range(len(histogram))]

    height, width = np.shape(image)
    for i in range(width):
        for j in range(height):
            x_components[image[j][i]].append(i)
            y_components[image[j][i]].append(j)

    return histogram, bin_edges, x_components, y_components


def find_marker(image: np.ndarray, histogram, x_components, y_components):
    print('Finding marker')
    ratio = []
    w, h = np.shape(image)

    for i in range(len(x_components)):
        if histogram[i] > 0:
            x_min = np.min(x_components[i])
            x_max = np.max(x_components[i])
            y_min = np.min(y_components[i])
            y_max = np.max(y_components[i])
            non_z = len(x_components[i])
            all_px = (x_max - x_min) * (y_max - y_min)
            square_to_circle_ratio = float(non_z) / float(all_px)

            # compensates for small noisy components and small excised roots
            if float(non_z) / float(w * h) > 0.0001:
                # the inscribed circle of a bounding box fills exactly 78.64 percent.
                # We allow 8.64 percent variation due to noise
                if square_to_circle_ratio > 0.7:
                    # sanity check
                    if (float(y_max) - float(y_min)) > 0:
                        # determine tag ratio
                        tagRatio = (float(x_max) - float(x_min)) / (float(y_max) - float(y_min))
                        ratio.append(np.abs(1 - (float(x_max) - float(x_min)) / (float(y_max) - float(y_min))))
                    else:
                        ratio.append(1000)
                else:
                    ratio.append(1000)
            else:
                ratio.append(1000)
        else:
            ratio.append(1000)

    rect = np.min(ratio)
    rect_idx = list(ratio).index(rect)

    x_min = np.min(x_components[rect_idx])
    x_max = np.max(x_components[rect_idx])
    y_min = np.min(y_components[rect_idx])
    y_max = np.max(y_components[rect_idx])

    idx = np.where(image == rect_idx)

    # bounding box
    i_min = np.min(idx[0])
    j_min = np.min(idx[1])
    i_max = np.max(idx[0])
    j_max = np.max(idx[1])

    sel = image != rect_idx
    image[sel] = 0

    if rect > 0.2:
        print('No marker detected')
        rect = 1
        rect_idx = 0

    return rect_idx, rect, i_max, i_min, j_min, j_max, image[i_min:i_max, j_min:j_max]


def find_crown(image: np.ndarray, histogram):
    print('Finding crown')
    image_copy = image.copy()
    height, width = np.shape(image_copy)
    found = False
    idx_1 = 0
    count = 0

    # We keep this piece of debug code, because the problem occurs in 1 of 10,000 images.
    # Perhaps we understand it one day.
    while not found:
        idx_1 = np.argmax(histogram)
        idx = np.where(image_copy == idx_1)
        if (np.max(idx[0]) + 1) == width and (np.max(idx[1]) + 1) == height and (np.min(idx[0])) == 0 and (
                np.min(idx[1])) == 0:
            if count < len(histogram):
                found = False
                count += 1
            else:
                found = True
            print('Only 1 background component that is smaller than the foreground ??? Probably a bug in the Masking routine')
        else:
            found = True

    # bounding box
    i_min = np.min(idx[0])
    i_max = np.max(idx[0])
    j_min = np.min(idx[1])
    j_max = np.max(idx[1])

    return idx_1, idx, i_max, i_min, j_min, j_max


# def find_tag(
#         image,
#         labelled_image,
#         masked_image,
#         circle_width,
#         circle_height,
#         histogram,
#         comps_x,
#         comps_y):
#     print('Finding tag')
#     ratio = []
#
#     for i in range(len(comps_x)):
#         if histogram[i] > 0:
#             x_min = np.min(comps_x[i])
#             x_max = np.max(comps_x[i])
#             y_min = np.min(comps_y[i])
#             y_max = np.max(comps_y[i])
#
#             # The tag should cover at least 0.5% of the picture
#             if 0.005 < float((x_max - x_min) * (y_max - y_min)) / float(circle_height * circle_width) < 0.01:
#                 # The tag should have more length then height
#                 if (x_max - x_min) >= 1.5 * (y_max - y_min):
#                     # The tag should be in the upper half of the image
#                     if y_min < (circle_height * 0.5):
#                         tagRatio = (float(x_max) - float(x_min)) / (float(y_max) - float(y_min))
#                         print('TagRatio detected: ' + str(tagRatio))
#                         ratio.append((float(x_max) - float(x_min)) / (float(y_max) - float(y_min)))
#                     else:
#                         ratio.append(-1)
#                 else:
#                     ratio.append(-1)
#             else:
#                 ratio.append(-1)
#         else:
#             ratio.append(-1)
#
#     rect = np.max(ratio)
#     if rect >= 0:
#         rectIdx = list(ratio).index(rect)
#         x_min = np.min(comps_x[rectIdx])
#         x_max = np.max(comps_x[rectIdx])
#         y_min = np.min(comps_y[rectIdx])
#         y_max = np.max(comps_y[rectIdx])
#     else:
#         x_min = x_max = y_min = y_max = 0
#         rectIdx = -1
#
#     print('Tag Ratio: ' + str(rect))
#     if rect == -1:
#         rectIdx = -1
#         i_min = i_max = j_min = j_max = 0
#     else:
#         idx = np.where(labelled_image == rectIdx)
#
#         # bounding box
#         i_min = np.min(idx[0])
#         j_min = np.min(idx[1])
#         i_max = np.max(idx[0])
#         j_max = np.max(idx[1])
#
#     sel = labelled_image != rectIdx
#     tag_crop = 10
#     if rect >= 0:
#         labelled_image[sel] = 0
#         try:
#             print('Checking for text')
#             tagText = ocr.find_text(
#                 image[i_min + tag_crop:i_max - tag_crop, j_min + tag_crop:j_max - tag_crop])
#         except:
#             tagText = 'Tag text extraction failed'
#             pass
#
#         # this order prefers the barcode over the text reader
#         try:
#             print('Check for barcode')
#             tagCode = ocr.find_barcode(
#                 image[i_min + tag_crop:i_max - tag_crop, j_min + tag_crop:j_max - tag_crop])
#             if len(tagCode) > 2:
#                 tagText = tagCode[8:len(tagCode) - 1]
#                 print(tagText)
#             else:
#                 print('Barcode too short: ' + tagCode)
#         except:
#             pass
#     else:
#         tagText = 'No label found'
#
#     return rectIdx, rect, i_max, i_min, j_min, j_max, masked_image[i_min + tag_crop:i_max - tag_crop, j_min + tag_crop:j_max - tag_crop], tagText


def find_tag(
        image,
        histogram,
        exclude_idx,
        root_idx,
        root_idx_list):
    print('Finding tag')
    image_copy = image.copy()
    counter = 0
    again = True

    while again:
        counter += 1
        if counter > 30:
            break

        again = False
        objects = len(histogram)
        if objects > 1:
            comp = np.argmax(histogram)
            excluded = 0
            if comp in exclude_idx:
                histogram[comp] = 0
                comp = np.argmax(histogram)
                if excluded > len(exclude_idx):
                    print('Error: Image is not usable')
                    break
                else:
                    excluded += 1
            idx2 = comp
        else:
            idx2 = -1

        image_copy = image.copy()
        idx = np.where(image_copy == idx2)

        # bounding box
        try:
            i_min = np.min(idx[0])
            i_max = np.max(idx[0])
            j_min = np.min(idx[1])
            j_max = np.max(idx[1])

            sel = image_copy != idx2
            image_copy[sel] = 0
            non_z = len(idx[0])
            bounding_box_size = (i_max - i_min) * (j_max - j_min)
            zeros = bounding_box_size - non_z
            ratio = float(zeros) / float(non_z)

            if counter >= objects:
                again = False
        except:
            again = True
    else:
        histogram[root_idx] = 0
        result = np.zeros_like(image_copy)
        result[root_idx_list] = 1
        return idx2, histogram, i_max, i_min, j_min, j_max  # [right:left, bottom:top], histogram


def find_excised(
        image: np.ndarray,
        histogram,
        exclude_idx):
    print('Finding excised root')
    image_copy = image.copy()
    counter = 0
    again = True

    while again:
        counter += 1
        if counter > 30:
            break

        again = False
        objects = len(histogram)
        if objects > 1:
            component = np.argmax(histogram)
            excluded = 0
            if component in exclude_idx:
                histogram[component] = 0
                component = np.argmax(histogram)
                if excluded > len(exclude_idx):
                    print("Unable to find excised root")
                    break
                else:
                    excluded += 1
            idx2 = component
        else:
            idx2 = -1

        image_copy = image.copy()
        idx = np.where(image_copy == idx2)

        # bounding box
        try:
            iMin = np.min(idx[0])
            jMin = np.min(idx[1])
            iMax = np.max(idx[0])
            jMax = np.max(idx[1])

            selection = image_copy != idx2
            image_copy[selection] = 0
            nonZ = len(idx[0])
            boundingBoxSize = (iMax - iMin) * (jMax - jMin)
            zeros = boundingBoxSize - nonZ
            ratio = float(zeros) / float(nonZ)

            if counter >= objects:
                again = False
        except:
            again = True

    selection = image_copy == idx2
    image_copy[selection] = 255

    return idx2, iMin, iMax, jMin, jMax, image_copy[iMin:iMax, jMin:jMax], (iMax + iMin) / 2, (jMax + jMin) / 2,


def medial_axis_skeleton(image: np.ndarray) -> (np.ndarray, np.ndarray):
    # to achieve consistent result between distance field and medial axis skeleton, set image borders to black
    image[0, :] = 0
    image[len(image) - 1, :] = 0
    image[:, len(image[0]) - 1] = 0
    image[:, 0] = 0

    distance = np.sqrt(distance_transform_edt(image > 0)) * 2
    med_axis = medial_axis(image > 0)

    return med_axis, distance


def label(image: np.ndarray) -> np.ndarray:
    labeled, nr_objects = ndimage.label(image)

    if nr_objects == 0:
        return None

    val = labeled.flatten()
    hist = []
    hist += range(np.max(val) + 1)
    test, _ = np.histogram(val, hist)
    comp1 = np.max(test)
    idx1 = list(test).index(comp1)

    if nr_objects > 1:
        test[idx1] = 0
        comp2 = np.max(test)
        idx2 = list(test).index(comp2)
        test[idx2] = 0
    else:
        idx2 = 1

    idx = np.where(labeled == idx2)

    # bounding box
    i_min = np.min(idx[0])
    j_min = np.min(idx[1])
    i_max = np.max(idx[0])
    j_max = np.max(idx[1])

    # just return the cropped image of the largest component
    return labeled[i_min:i_max, j_min:j_max]


def segment_internal(options: DIRTOptions, image: np.ndarray):
    output_prefix = join(options.output_directory, options.input_stem)

    # find features
    labeled, _ = ndimage.label(image)

    # calculate histogram
    histogram, bin_edges, x_components, y_components = calculate_histogram(labeled)

    # find crown
    crown, crown_list, crown_left, crown_right, crown_bottom, crown_top = find_crown(
        labeled.copy(),
        histogram)

    # find marker
    marker, marker_ratio, marker_left, marker_right, marker_top, marker_bottom, _ = find_marker(
        labeled.copy(),
        histogram,
        x_components,
        y_components)
    marker_width = float(marker_left) - float(marker_right)
    marker_height = float(marker_bottom) - float(marker_top)

    try:
        histogram[marker] = 0
    except:
        pass

    # find tag
    tag, tag_list, tag_left, tag_right, tag_bottom, tag_top = find_tag(
        labeled.copy(),
        histogram,
        [marker, crown],
        crown,
        crown_list)

    try:
        histogram[tag] = 0
    except:
        pass

    # find excised root(s)
    # if options.excised_roots > 1 and options.crown_root:
    #     for i in range(options.excised_roots):
    #         exRIdx, excised_left, excised_right, excised_top, excised_bottom, imgExRoot, centerPtx, centerPty = find_excised(
    #             labeled.copy(),
    #             histogram,
    #             [marker, tag, crown])
    #         print(f"Found excised root {i}")
    # elif options.excised_roots == 1 and not options.crown_root:
    #     exRIdx, excised_left, excised_right, excised_top, excised_bottom, imgExRoot, centerPtx, centerPty = find_excised(
    #         labeled.copy(),
    #         histogram,
    #         [marker, tag])

    # save masked image to file
    imageio.imwrite(f"{output_prefix}.mask.png", skimage.img_as_uint(image))

    # save crown bounding box overlay to file
    fig, ax = plt.subplots(1)
    ax.imshow(image)
    crown_patch = patches.Rectangle(
        (crown_top, crown_left),
        crown_bottom - crown_top,
        crown_right - crown_left,
        edgecolor='r')
    crown_patch.set_alpha(0.5)
    ax.add_patch(crown_patch)
    plt.axis('off')
    plt.savefig(f"{output_prefix}.bounding.crown.png")
    plt.clf()

    # save marker bounding box overlay to file
    fig, ax = plt.subplots(1)
    ax.imshow(image)
    marker_patch = patches.Rectangle(
        (marker_top, marker_left),
        marker_bottom - marker_top,
        marker_right - marker_left,
        edgecolor='r')
    marker_patch.set_alpha(0.5)
    ax.add_patch(marker_patch)
    plt.axis('off')
    plt.savefig(f"{output_prefix}.bounding.marker.png")
    plt.clf()

    # save tag bounding box overlay to file
    fig, ax = plt.subplots(1)
    ax.imshow(image)
    tag_patch = patches.Rectangle(
        (tag_top, tag_left),
        tag_bottom - tag_top,
        tag_right - tag_left,
        edgecolor='r')
    tag_patch.set_alpha(0.5)
    ax.add_patch(tag_patch)
    plt.axis('off')
    plt.savefig(f"{output_prefix}.bounding.tag.png")
    plt.clf()

    return '', marker_ratio, marker_width, marker_height


def segment(options: DIRTOptions, image: np.ndarray) -> DIRTResults:
    print(f"Segmenting file: {options.input_name}")
    start = time.time()

    tag, marker_ratio, marker_width, marker_height = segment_internal(options=options, image=image)
    x_scale = options.marker_diameter / float(marker_width)
    y_scale = options.marker_diameter / float(marker_height)

    seconds = time.time() - start
    print(f"Segmentation finished in {seconds} seconds")

    return DIRTResults(
        image_name=options.input_stem,
        tag=tag,
        marker_ratio=round(float(marker_ratio), 3),
        marker_width=marker_width,
        marker_height=marker_height,
        x_scale=1.0 if x_scale <= 0.0 else x_scale,
        y_scale=1.0 if y_scale <= 0.0 else y_scale)
