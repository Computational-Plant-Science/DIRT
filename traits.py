import time
from math import ceil, sqrt
from os.path import join
from typing import List, Literal

import imageio
import numpy as np
import warnings
import scipy.stats
import scipy.misc
import scipy.interpolate
import matplotlib.pyplot as plt
from matplotlib import patches
from skimage import img_as_uint
import graph_tool.draw as gt
import graph_tool.centrality as gc
import graph_tool.topology as gtop
from graph_tool import Graph, GraphView
from graph_tool.libgraph_tool_core import VertexBase

import kmeans as km
import ransac as rs
from skeletonization import extract_graph, initial_branching_point, thickest_full_path, root_tip_skeleton, \
    hypocotol_cluster, cluster_overlay_image, lateral_lengths
from options import DIRTOptions
from results import DIRTResults
from segmentation import medial_axis_skeleton

warnings.simplefilter('ignore', np.RankWarning)


def filter_root_tip_path_tangent(path: List[VertexBase], branch_points, window=5):
    tangents = []
    for point in branch_points:
        tangent = []
        try:
            idx = path.index(point)
            for j in range(window):
                if j < idx:
                    try:
                        tangent.append(path[idx - j])
                    except:
                        pass
                if j < len(path) - j:
                    try:
                        tangent.append(path[idx + j])
                    except:
                        pass
        except:
            pass
        tangents.append(tangent)

    return tangents


def roots_per_segment(adventitious_count, basal_count, cDiaAdv, cDiaBas):
    adv_roots = 0
    bas_roots = 0
    seg1_roots = adventitious_count[0] - adventitious_count[len(adventitious_count) - 1]
    seg2_roots = basal_count[0] - basal_count[len(basal_count) - 1]

    for idx, i in enumerate(adventitious_count):
        if adventitious_count[idx - 1] - i > 0:  # just a little threshold to get rid of noise
            adv_roots += 1
    for idx, i in enumerate(basal_count):
        if basal_count[idx - 1] - i > 0:  # just a little threshold to get rid of noise
            bas_roots += 1

    return adv_roots, bas_roots, seg1_roots, seg2_roots, np.mean(cDiaAdv), np.mean(cDiaBas)


def fit_line(graph: Graph, path: List[VertexBase]):
    print("Fitting polynomial")
    if len(path) == 0:
        return -1, -1

    x = []
    y = []
    vp = graph.vertex_properties["vp"]
    for vertex in path:
        x.append(vp[graph.vertex(vertex)]['coord'][0])
        y.append(vp[graph.vertex(vertex)]['coord'][1])

    return np.polyfit(x, y, 1)


def fit_line_xy(x, y, ransac=False):
    print("Fitting polynomial")
    if ransac:
        try:
            x, y = rs.ransac_fit(x, y)
        except:
            print("RANSAC fitting failed, using simple linear fitting")
    try:
        return np.polyfit(x, y, 1)
    except:
        print("Failed to fit polynomial")
        return 1.0, 1.0


def angle_to_x_ax(graph: Graph, path: List[VertexBase]):
    m2, _ = fit_line(graph, path)
    alpha = (np.arctan(m2) * 180) / np.pi
    return np.fabs(alpha)


def angle_to_x_ax_xy(x, y, ransac=False):
    m2, _ = fit_line_xy(x, y, ransac)
    alpha = (np.arctan(m2) * 180) / np.pi
    return np.fabs(alpha)


def angle_at_distance(
        graph: Graph,
        path: List[VertexBase],
        lat,
        branch_points,
        scale: float,
        distance=20):
    """
    Calculates the angle between the hypotocol and all root tip paths.
    """

    angles = []
    scaled_distance = int(distance / scale)
    for i in range(len(lat)):
        if len(lat[i]) > scaled_distance and graph.vertex(branch_points[i]) in path:
            angle = angle_to_x_ax_xy(graph, lat[i][:scaled_distance])
            angles.append(angle)

    return np.mean(angles)


def angles_per_cluster(
        graph: Graph,
        path: List[VertexBase],
        cAdv,
        cBas,
        lat,
        branch_points,
        scale,
        dist=20):
    """
    Estimate angles per cluster to distinguish adventitious roots from basal roots.
    """

    vp = graph.vertex_properties["vp"]
    minPathsAdv = np.min(cAdv)
    for idx, i in enumerate(path):
        if minPathsAdv > vp[i]['nrOfPaths']:
            mean_adventitious = angle_at_distance(graph, path[:idx], lat, branch_points, scale, dist)
            mean_basal = angle_at_distance(graph, path[idx:], lat, branch_points, scale, dist)
            break

    return mean_adventitious, mean_basal


def path_length(path: List[VertexBase], scale: float):
    return len(path) * scale


def angle_between_paths(graph: Graph, path1: List[VertexBase], path2: List[VertexBase]):
    """
    Calculates the angle between two paths.
    """

    m1, _ = fit_line(graph, path1)
    m2, _ = fit_line(graph, path2)
    up = m1 - m2
    low = 1 + (m1 * m2)
    tan_alpha = np.fabs(up / low)
    alpha = (np.arctan(tan_alpha) * 180) / np.pi

    return np.fabs(alpha)


def lateral_angles(
        graph: Graph,
        path: List[VertexBase],
        lat,
        corrBranchpts):
    """
    Calculates angles between the hypocotol and all paths
    """

    angles = []
    tangents = filter_root_tip_path_tangent(path, corrBranchpts)
    for i in range(len(lat)):
        angle = angle_between_paths(graph, tangents[i], lat[i])
        angles.append(angle)
    try:
        min_angle = np.min(angles)
        max_angle = np.max(angles)
        angle_range = max_angle - min_angle
        angle_avg = np.average(angles)
    except:
        min_angle = -1
        max_angle = -1
        angle_range = -1
        angle_avg = -1

    return angle_range, angle_avg, min_angle, max_angle, angles


def symmetry(graph: Graph, root_tips: List[VertexBase]):
    x_max = 0
    x_min = 5000000000
    y_max = 0
    y_min = 5000000000
    x_sum = 0
    y_sum = 0
    count = 0.0
    vp = graph.vertex_properties["vp"]

    # calculate bounding box
    for tip in root_tips:
        # for i in r:
        # v = graph.vertex(i)
        if vp[tip]['coord'][0] > x_max:
            x_max = vp[tip]['coord'][0]
        if vp[tip]['coord'][0] < x_min:
            x_min = vp[tip]['coord'][0]
        if vp[tip]['coord'][1] > y_max:
            y_max = vp[tip]['coord'][1]
        if vp[tip]['coord'][1] < y_min:
            y_min = vp[tip]['coord'][1]
        x_sum += vp[tip]['coord'][0]
        y_sum += vp[tip]['coord'][1]
        count += 1.0

    if count > 0:
        x_avg = x_sum / count
        y_avg = y_sum / count
        x_avg_box = (x_max - x_min) / 2
        y_avg_box = (y_min - y_max) / 2
        x_vec = x_avg_box - x_avg
        y_vec = y_avg_box - y_avg
        symmetry = [x_vec, y_vec]
    else:
        symmetry = [np.nan, np.nan]

    return symmetry


def mean_lateral_length(self, graph: Graph, paths: List[VertexBase]):
    x = []
    lengths = []
    vp = graph.vertex_properties["vp"]

    if len(paths) > 0:
        for i in paths:
            length = len(i)
            if len(i) > 0:
                lengths.append(length)
                x.append(vp[graph.vertex(i[0])]['imgIdx'][0])
            else:
                lengths.append(-1)
                x.append(-1)

    f2 = scipy.interpolate.interp1d(x, lengths, kind='cubic')
    return lengths, x, f2(x), np.average(lengths) * self.__scale


def path_diameter_quantiles(graph, path, scale: float):
    filtered_graph = filter_path_diameters(graph.copy(), path, scale)
    x = []
    y = []
    length = 0
    vp = filtered_graph.vertex_properties["vp"]
    for i in path:
        length += 1
        if vp[i]['diameter'] > 0:
            x.append(length)
            y.append(vp[i]['diameter'])

    l = len(y) - 1
    l25 = int(l * 0.25)
    l50 = int(l * 0.5)
    l75 = int(l * 0.75)
    l90 = int(l * 0.90)

    d25 = np.average(y[:l25])
    d50 = np.average(y[l25:l50])
    d75 = np.average(y[l50:l75])
    d90 = np.average(y[l90:])

    return d25, d50, d75, d90, x, y


def filter_path_diameters(graph: Graph, path: List[VertexBase], scale: float):
    """
    Remove diameters around branching points in range of the branching point diameter.
    """

    vp = graph.vertex_properties["vp"]
    for i in range(len(path)):
        count = 0
        for _ in path[i].out_neighbours():
            count += 1
            if count > 2:
                break
        if count > 2:
            for j in range(int(vp[path[i]]['diameter'] / scale)):
                if j > 20: break
                if i - j > 0: vp[path[i - j]]['diameter'] = 0
                if i + j < len(path): vp[path[i + j]]['diameter'] = 0
            vp[path[i]]['diameter'] = 0

    return graph


def smooth(
        x,
        window_len=11,
        window_type: Literal['flat', 'hanning', 'hamming', 'bartlett', 'blackman'] = 'hanning'):
    """Smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise (ValueError, "Only 1-dimension arrays are allowed.")

    if x.size < window_len:
        raise (ValueError, "Array must be larger than window size.")

    if window_len < 3:
        return x

    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    w = np.ones(window_len, 'd') if window_type == 'flat' else eval('np.' + window_type + '(window_len)')
    y = np.convolve(w / w.sum(), s, mode='valid')

    return y


def histogram_peaks(angle):
    print('Finding root angle histogram peaks')
    try:
        pdf, _ = np.histogram(angle, bins=9, range=(0, 90))
        maximum = []
        for i in range(0, len(pdf) - 1):
            if pdf[i - 1] <= pdf[i] or i == 0:
                if pdf[i + 1] <= pdf[i] or i == len(pdf) - 1:
                    maximum.append(i)

        sorted_pdf = np.argsort(pdf)
        strongest_max = []
        for i in sorted_pdf[::-1]:
            if i in maximum:
                strongest_max.append(i)
            if len(strongest_max) == 2:
                break

        avg_result = []
        for i in strongest_max:
            avg_angle = []
            for j in angle:
                if j >= i * 10:
                    if j <= (i + 1) * 10:
                        avg_angle.append(j)
            avg_result.append(np.mean(avg_angle))

        if len(avg_result) > 1:
            return np.max(avg_result), np.min(avg_result)
        else:
            return avg_result[0], -1
    except:
        print('No histogram peaks found')
        print(pdf)
        return -1, -1


def paths_per_depth(graph: Graph, path: List[VertexBase]):
    print('Calculating root paths per depth')
    depth = [len(path)]
    paths = [len(path)]
    vp = graph.vertex_properties['vp']
    paths_50 = vp[path[0]]['nrOfPaths'] / 2
    drop_50 = 0
    for i in path:
        paths.append(vp[i]['nrOfPaths'])
        depth.append(vp[i]['coord'][1])
        if vp[i]['nrOfPaths'] <= paths_50:
            drop_50 = vp[i]['coord'][1]

    return float(drop_50), paths, depth


def diameter_radii(paths, dia, clusters):
    pts = []
    for i in range(len(paths)):
        pts.append([paths[i], dia[i]])

    cl = km.KMeans(pts)
    c = cl.kmeans(clusters, 0.01)

    cX = []
    cY = []

    for nc in range(clusters):
        cX.append([])
        cY.append([])
    for nc in range(clusters):
        for i in c[nc]:
            cX[nc].append(i[0])
            cY[nc].append(i[1])

    if clusters == 2:
        if cY[0][0] > cY[1][0]:
            return cX[0], cY[0], cX[1], cY[1]
        else:
            return cX[1], cY[1], cX[0], cY[0]
    if clusters == 3:
        if cY[0][0] > cY[1][0] > cY[2][0]:
            return cX[0], cY[0], cX[1], cY[1], cX[2], cY[2]
        if cY[0][0] > cY[1][0] and cY[1][0] < cY[2][0]:
            return cX[0], cY[0], cX[2], cY[2], cX[1], cY[1]
        if cY[0][0] < cY[1][0] < cY[2][0]:
            return cX[2], cY[2], cX[1], cY[1], cX[0], cY[0]
        else:
            return cX[2], cY[2], cX[0], cY[0], cX[1], cY[1]


def root_top_angle(graph: Graph, lat):
    """
    Calculate angle between the hypocotyl and all root tip paths
    """

    angles = []
    for i in range(len(lat)):
        ang = angle_to_x_ax(graph, lat[i])
        angles.append(ang)

    if len(angles) != 0:
        min = np.min(angles)
        max = np.max(angles)
        rng = max - min
        med = np.median(angles)
        return med, min, max, rng, angles
    else:
        return -1, -1, -1, -1, angles


def filter_cdf_tangent_slope(x, cdf, idx, window=20):
    """
    Estimate the slope at a point over a small region along the CDF.
    """

    tangents = []
    for i in idx:
        tangent_x = []
        tangent_y = []
        for j in range(window):
            if j < i:
                try:
                    tangent_x.append(x[int(i) - j])
                    tangent_y.append(cdf[int(i) - j])
                except:
                    pass
            if j < len(cdf) - j:
                try:
                    tangent_x.append(x[int(i) + j])
                    tangent_y.append(cdf[int(i) + j])
                except:
                    pass
        a, _ = fit_line_xy(tangent_x, tangent_y)
        tangents.append(a)

    return tangents


def mask_traits(image: np.ndarray, x_scale: float, y_scale: float):
    """
    Compute all mask-based traits at once.
    """

    # We compute here all mask based traits at once (long lives Spagetti code :-) )
    if len(image) > 0:
        print('IMG OK')
    else:
        print('IMAGE IS NOT ACCESSIBLE')
        return 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', ['nan'] * 9, ['nan'] * 9, 'nan'
    h, w = np.shape(image)
    xx = []
    yy = []
    count = 0
    black_pixels = []
    root_densities = []

    # compute density value
    for i in range(h):
        idx = np.where(image[i] > 0)
        try:
            start = idx[0][0]
            white = len(idx[0])
            end = idx[0][len(idx[0]) - 1]
            black = (end - start) - white
            count += white
            width = float(end - start)
            xx.append(float(i) * y_scale)
            if black == 0: black = 1
            black_pixels.append(black)
            yy.append(float(width) * x_scale)

            if w > 0:
                normWhite = float(white) / float(w)
                if normWhite == 1.: normWhite = 0.
                if black > 0:
                    normBlack = float(black) / float(w)
                else:
                    normBlack = 1.

                if normWhite > 0.: root_densities.append(float(normWhite) / float(normBlack))
        except:
            root_densities.append(float(0.0))
            print('empty image line in crown file -> placed 0. as density for this line')
            pass

    avg_root_density = np.average(root_densities)
    print('Avg. Root density: ' + str(avg_root_density))

    y_smooth = yy
    smooth_region = 15
    xx_norm = np.array(xx) / np.max(xx)
    ten_percent = float(len(yy)) * 0.1

    # retrieve stem diameter as the average of the distance field in the first 10% of height
    try:
        stem_diameter = np.median(xx[20:int(ten_percent)])
    except:
        stem_diameter = -1

    # compute a simple angle at top and bottom along the outline for monocots.
    # this is noisier than the D10 or D20 values, which are more robust.
    try:
        top_angle = angle_to_x_ax_xy(
            xx[int(ten_percent):int(ten_percent) * 3],
            yy[int(ten_percent):int(ten_percent) * 3],
            True)
        bottom_angle = angle_to_x_ax_xy(
            xx[int(ten_percent) * 3:int(ten_percent) * 9],
            yy[int(ten_percent) * 3:int(ten_percent) * 9],
            True)
    except:
        top_angle = -1
        bottom_angle = -1

    # smooth the noisy data
    for i in range(smooth_region, len(yy) - smooth_region):
        tmp = []
        for j in range(smooth_region):
            tmp.append(yy[i + j])
            tmp.append(yy[i - j])
        y_smooth[i] = np.median(tmp)

    # compute the width parameters of the root
    width_median = np.median(y_smooth)
    width_max = np.max(y_smooth)
    # compute the cummulative width profile
    y_cwp = np.array(y_smooth).cumsum()
    y_cwp = y_cwp / np.max(y_cwp)

    # find index at x% width accumulation
    D = []
    dD = 0.1
    count = 0.0
    for i in y_cwp:
        count += 1.0
        if i > dD:
            D.append(count)
            dD += 0.1
        if dD == 1.0:
            break
    while len(D) < 9:
        D.append(-1)

    d_slope = filter_cdf_tangent_slope(xx_norm, y_cwp, D[:len(D) - 1], 100)
    # convert D array counts to percentages saved in the output file
    D = np.array(D, dtype=float) / float(len(y_cwp))

    return DIRTResults(
        mean_density=round(float(avg_root_density), 3),
        median_width=round(float(width_median), 3),
        max_width=round(float(width_max), 3),
        acc_width_10=round(float(D[0]), 3),
        acc_width_20=round(float(D[1]), 3),
        acc_width_30=round(float(D[2]), 3),
        acc_width_40=round(float(D[3]), 3),
        acc_width_50=round(float(D[4]), 3),
        acc_width_60=round(float(D[5]), 3),
        acc_width_70=round(float(D[6]), 3),
        acc_width_80=round(float(D[7]), 3),
        acc_width_90=round(float(D[8]), 3),
        acc_width_slope_10=round(float(d_slope[0]), 3),
        acc_width_slope_20=round(float(d_slope[1]), 3),
        acc_width_slope_30=round(float(d_slope[2]), 3),
        acc_width_slope_40=round(float(d_slope[3]), 3),
        acc_width_slope_50=round(float(d_slope[4]), 3),
        acc_width_slope_60=round(float(d_slope[5]), 3),
        acc_width_slope_70=round(float(d_slope[6]), 3),
        acc_width_slope_80=round(float(d_slope[7]), 3),
        acc_width_slope_90=round(float(d_slope[8]), 3),
        max_diameter_90=count * (x_scale * y_scale),
        # drop_50=stem_diameter,
        top_angle=round(float(top_angle), 3),
        bottom_angle=round(float(bottom_angle), 3))


def angle_quantiles(graph: Graph, lat):
    angles_25 = []
    angles_50 = []
    angles_75 = []
    angles_90 = []

    for i in range(len(lat)):
        length = len(lat[i])
        length_25 = int(length * 0.25)
        length_50 = int(length * 0.5)
        length_75 = int(length * 0.75)
        length_90 = int(length * 0.90)
        angle_25 = angle_to_x_ax(graph, lat[i][:length_25])
        angle_50 = angle_to_x_ax(graph, lat[i][:length_50])
        angle_75 = angle_to_x_ax(graph, lat[i][:length_75])
        angle_90 = angle_to_x_ax(graph, lat[i][:length_90])
        angles_25.append(angle_25)
        angles_50.append(angle_50)
        angles_75.append(angle_75)
        angles_90.append(angle_90)

    return angles_25, angles_50, angles_75, angles_90


def path_diameter(graph: Graph, path: List[VertexBase]):
    vp = graph.vertex_properties["vp"]
    x = []
    y = []
    length = 0

    for i in path:
        length += 1
        if vp[i]['diameter'] > 0:
            x.append(length)
            y.append(vp[i]['diameter'])

    coeffs = np.polyfit(x, y, 1)
    mean_diameter = np.average(y)

    return mean_diameter, coeffs[0]


def path_branching_frequency(
        root_tips: List[VertexBase],
        path: List[VertexBase],
        scale: float):
    branch_points = []
    for tip in root_tips:
        if tip[0] != path[0]:
            branch_points.append(tip[0])

    try:
        branching_frequency = float(len(path)) / float(len(np.unique(branch_points)))
    except:
        branching_frequency = -1

    return branching_frequency * scale


def traits(options: DIRTOptions, image: np.ndarray, results: DIRTResults) -> DIRTResults:
    start_traits = time.time()
    output_prefix = join(options.output_directory, options.input_stem)

    # mask traits
    start = time.time()
    x_scale = results['x_scale']
    y_scale = results['y_scale']
    results = {**results, **mask_traits(image, x_scale, y_scale)}
    print(f"Mask traits computed in {ceil(time.time() - start)} seconds")

    # medial axis
    start = time.time()
    med_axis, med_axis_diameter = medial_axis_skeleton(image)
    imageio.imwrite(f"{output_prefix}.medial.png", img_as_uint(med_axis))
    med_axis_img = img_as_uint(med_axis)
    med_axis_overlay = image.copy()
    med_axis_overlay[np.where(med_axis_img != 0)] = 120
    imageio.imwrite(f"{output_prefix}.medial.overlay.png", med_axis_overlay)
    print(f"Medial axis computed in {ceil(time.time() - start)} seconds")

    # thickest path
    start = time.time()
    graph = extract_graph(image, med_axis_diameter, x_scale, y_scale)
    num_branching_pts = graph.num_vertices()
    top_branch_point = initial_branching_point(graph, np.shape(med_axis)[0])
    thickest_path = thickest_full_path(graph, top_branch_point)
    vp = graph.vertex_properties["vp"]
    thickest_path_overlay = image.copy()
    thickest_path_overlay[np.where(image != 0)] = 120
    for point in thickest_path:
        thickest_path_overlay[int(vp[point]['abs_coord'][1]), int(vp[point]['abs_coord'][0])] = 255
    imageio.imwrite(f"{output_prefix}.thickest.overlay.png", thickest_path_overlay)
    current = DIRTResults(branching_points=num_branching_pts)
    results = {**results, **current}
    print(f"Thickest path computed in {ceil(time.time() - start)} seconds")

    # minimum spanning tree
    # start = time.time()
    # mst = GraphView(graph, efilt=gtop.min_spanning_tree(graph))
    # pos = gt.sfdp_layout(mst)
    # gt.graph_draw(
    #     mst,
    #     pos=pos,
    #     edge_pen_width=mst.edge_properties['w'],
    #     output=f"{output_prefix}.mst.png")
    # print(f"Minimum spanning tree computed in {ceil(time.time() - start)} seconds")

    # skeleton
    start = time.time()
    skeleton, paths, tips, tip_median_diameter, tip_mean_diameter, \
    skeleton_depth, skeleton_width, skeleton_height = root_tip_skeleton(graph, thickest_path)
    current = DIRTResults(
        tip_median_diameter=round(float(max(tip_median_diameter)), 3),
        tip_mean_diameter=round(float(tip_mean_diameter), 3),
        paths=len(paths),
        depth=round(float(skeleton_depth), 3),
        width=round(float(skeleton_width), 3))
    results = {**results, **current}
    print(f"Skeleton, paths, and tips computed in {ceil(time.time() - start)} seconds")

    # spatial distribution
    start = time.time()
    x_distribution, y_distribution = symmetry(skeleton, tips)
    current = DIRTResults(spatial_dist_x=x_distribution, spatial_dist_y=y_distribution)
    results = {**results, **current}
    print(f"Spatial distribution computed in {ceil(time.time() - start)} seconds")

    # hypocotyl
    start = time.time()
    branchRad, nrPaths = hypocotol_cluster(skeleton, thickest_path)
    print(f"Hypocotol computed in {ceil(time.time() - start)} seconds")

    # diameter radii
    start = time.time()
    c1x, c1y, c2x, c2y = diameter_radii(nrPaths, branchRad, 2)
    print(f"2 k-means clusters computed in {ceil(time.time() - start)} seconds")

    # feature classes
    start = time.time()
    overlay = cluster_overlay_image(skeleton, thickest_path, image, x_scale, y_scale, c1x, c1y, c2x, c2y)
    adventitious, basal, path_seg1, path_seg2, hypocotyl_diameter, taproot_diameter = roots_per_segment(c1y, c2y, c1x, c2x)
    current = DIRTResults(
        adventitious=adventitious,
        basal=basal,
        path_seg_1=path_seg1,
        path_seg_2=path_seg2,
        hypocotol_diameter=round(float(hypocotyl_diameter), 3),
        taproot_diameter=round(float(taproot_diameter), 3))
    results = {**results, **current}
    print(f"Root classes computed in {ceil(time.time() - start)} seconds")

    # lateral lengths
    start = time.time()
    laterals, corr_branch_pts = lateral_lengths(skeleton, tips, (x_scale + y_scale) / 2)
    print(f"Lateral lengths computed in {ceil(time.time() - start)} seconds")

    start = time.time()
    if c1x is not None and c1y is not None and c2x is not None and c2y is not None:
        adventitious_angles, basal_angles = angles_per_cluster(
            skeleton,
            thickest_path,
            c1y,
            c2y,
            laterals,
            corr_branch_pts,
            (x_scale + y_scale) / 2,
            dist=20)
    else:
        adventitious_angles = basal_angles = 'nan'
    current = DIRTResults(adventitious_angles=float(adventitious_angles), basal_angles=float(basal_angles))
    results = {**results, **current}
    print(f"Lateral angles at 2cm computed in {ceil(time.time() - start)} seconds")

    # angle quantiles
    start = time.time()
    try:
        a25, a50, a75, a90 = angle_quantiles(skeleton, laterals)
        print(f"Angle quantiles computed in {ceil(time.time() - start)} seconds")
    except:
        a25 = ['nan']
        a50 = ['nan']
        a75 = ['nan']
        a90 = ['nan']
        print('No angle quantile calculated')

    # top angle
    start = time.time()
    top_median, top_min, top_max, top_range, top_angles = root_top_angle(skeleton, laterals)
    current = DIRTResults(
        median_angle=float(top_median),
        min_angle=top_min,
        max_angle=top_max,
        lateral_angular_range=top_range)
    results = {**results, **current}
    print(f"Top angle characteristics computed in {ceil(time.time() - start)} seconds")

    # Soil tissue angle
    start = time.time()
    sta_range, sta_median, sta_min, sta_max, sta_angles = lateral_angles(
        graph,
        thickest_path,
        laterals,
        corr_branch_pts)
    current = DIRTResults(
        soil_tissue_angle_range=round(float(sta_range), 3),
        soil_tissue_median_angle=round(float(sta_median), 3),
        soil_tissue_min_angle=round(float(sta_min), 3),
        soil_tissue_max_angle=round(float(sta_max), 3))
    results = {**results, **current}
    print(f"Lateral angle characteristics computed in {ceil(time.time() - start)} seconds")

    # taproot quantiles
    start = time.time()
    tr_dia_25, tr_dia_50, tr_dia_75, tr_dia_90, dqx, dqy = path_diameter_quantiles(
        graph,
        thickest_path,
        (x_scale + y_scale) / 2)
    current = DIRTResults(
        taproot_diameter_25=round(float(tr_dia_25), 3),
        taproot_diameter_50=round(float(tr_dia_50), 3),
        taproot_diameter_75=round(float(tr_dia_75), 3),
        taproot_diameter_90=round(float(tr_dia_90), 3),
        soil_tissue_max_angle=round(float(sta_max), 3))
    results = {**results, **current}
    print(f"Diameter quantiles computed in {ceil(time.time() - start)} seconds")

    # dominant lateral angles
    start = time.time()
    sta_dom_1, sta_dom_2 = histogram_peaks(sta_angles)
    current = DIRTResults(
        soil_tissue_dominant_angle_1=round(float(sta_dom_1), 3),
        soil_tissue_dominant_angle_2=round(float(sta_dom_2), 3))
    results = {**results, **current}
    print(f"Dominant lateral angles computed in {ceil(time.time() - start)} seconds")

    # dominant lateral angle quantiles
    start = time.time()
    sta_dom_1_25, sta_dom_2_25 = histogram_peaks(a25)
    sta_dom_1_50, sta_dom_2_50 = histogram_peaks(a50)
    sta_dom_1_75, sta_dom_2_75 = histogram_peaks(a75)
    sta_dom_1_90, sta_dom_2_90 = histogram_peaks(a90)
    current = DIRTResults(
        soil_tissue_dominant_angle_1_25=round(float(sta_dom_1_25), 3),
        soil_tissue_dominant_angle_2_25=round(float(sta_dom_2_25), 3),
        soil_tissue_dominant_angle_1_50=round(float(sta_dom_1_50), 3),
        soil_tissue_dominant_angle_2_50=round(float(sta_dom_2_50), 3),
        soil_tissue_dominant_angle_1_75=round(float(sta_dom_1_75), 3),
        soil_tissue_dominant_angle_2_75=round(float(sta_dom_2_75), 3),
        soil_tissue_dominant_angle_1_90=round(float(sta_dom_1_90), 3),
        soil_tissue_dominant_angle_2_90=round(float(sta_dom_2_90), 3))
    results = {**results, **current}
    print(f"Dominant lateral angle quantiles computed in {ceil(time.time() - start)} seconds")

    # dominant top angles
    start = time.time()
    top_dom_1, top_dom_2 = histogram_peaks(top_angles)
    current = DIRTResults(top_dominant_angle_1=round(float(top_dom_1), 3), top_dominant_angle_2=round(float(top_dom_2), 3))
    results = {**results, **current}
    print(f"Dominant top angles computed in {ceil(time.time() - start)} seconds")

    print(f"All traits computed in {ceil(time.time() - start_traits)} seconds")
    return results
