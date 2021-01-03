import time
from collections import OrderedDict
from math import ceil

import imageio
import numpy as np
import skimage

from options import DIRTOptions
from results import DIRTResults
from traits import paths_per_depth, path_length, histogram_peaks, \
    path_diameter_quantiles, mean_lateral_length, angles_per_cluster, diameter_radii, \
    root_top_angle, lateral_angles, mask_traits, angle_quantiles, \
    path_diameter, path_branching_frequency
from skeletonization import extract_graph, \
    lateral, thickest_full_path, initial_branching_point, hypocotol_cluster, cluster_overlay_image, \
    lateral_lengths, root_tip_skeleton
from segmentation import label


# def process(
#         options: DIRTOptions,
#         results: DIRTResults) -> DIRTResults:
#
#     if options.excised_roots == 0:
#         return results
#
#     rtp_skel = -1
#     lateral_file_path = f"{options.input_stem}.lateral.png"
#     rtp = RootTips(options)
#     analyze = Analysis(options, (preprocessing_result.x_scale + preprocessing_result.y_scale) / 2)
#     image = imageio.imread(lateral_file_path, as_gray=True)
#     segmentation = Segmentation(image)
#     image_label = label(image)
#
#     if image_label is None:
#         raise ValueError(f"No image label found, skipping lateral file processing")
#     print(f"Processing lateral file: {lateral_file_path}")
#
#     skeleton, skeleton_diameter = find_medial_axis(image_label)
#     skeleton_graph = extract_graph(skeleton, skeleton_diameter, preprocessing_result.x_scale,
#                                    preprocessing_result.y_scale)
#     skeleton_root_vertex = lateral(skeleton_graph, np.shape(skeleton)[0])
#     thickest_path = thickest_path(skeleton, skeleton_root_vertex)
#
#     if options.traits.rtp_skeleton:
#         start = time.time()
#         rtp_skel, _, results['LT_MED_DIA'], results[
#             'LT_AVG_DIA'], _, rtps, _, _, _ = rtp.root_tip_skeleton(
#             thickest_path, skeleton_graph, True)
#         print(f"RTP skeleton computed in {ceil(time.time() - start)} seconds")
#
#     if rtp_skel != -1:
#         start = time.time()
#         results['LT_BRA_FRQ'] = path_branching_frequency(rtps, thickest_path)
#         results['NODAL_AVG_DIA'], _ = path_diameter(rtp_skel, thickest_path)
#         results['NODAL_LEN'] = path_length(
#             thickest_path,
#             (preprocessing_result.x_scale + preprocessing_result.y_scale) / 2)
#         print(f"Branching frequency computed in {ceil(time.time() - start)} seconds")
#
#     if rtp_skel != -1 and options.traits.rtp_lateral:
#         start = time.time()
#         lat, corrBranchpts, results['LT_DIST_FIRST'] = segmentation.findLaterals(
#             rtps,
#             rtp_skel,
#             (preprocessing_result.x_scale + preprocessing_result.y_scale) / 2,
#             thickest_path)
#         print(f"RTP laterals computed in {ceil(time.time() - start)} seconds")
#
#     if rtp_skel != -1 and options.traits.lateral_length:
#         start = time.time()
#         cmll_length, cmll_x, cmll_f2x, results['LT_AVG_LEN'] = mean_lateral_length(rtp_skel, lat, thickest_path)
#         np.savetxt(f"{options.input_stem}.length.histo.txt", cmll_length, delimiter=',')
#         np.savetxt(f"{options.input_stem}.length.x.txt", cmll_x, delimiter=',')
#         np.savetxt(f"{options.input_stem}.length.y.txt", cmll_f2x, delimiter=',')
#         print(f"Lateral length computed in {ceil(time.time() - start)} seconds")
#
#     if rtp_skel != -1 and options.traits.lateral_angles:
#         start = time.time()
#         results['LT_ANG_RANGE'], results['LT_AVG_ANG'], results['LT_MIN_ANG'], results[
#             'LT_MAX_ANG'], _ = lateral_angles(rtp_skel, thickest_path, lat, corrBranchpts)
#         print(f"Lateral angles computed in {ceil(time.time() - start)} seconds")