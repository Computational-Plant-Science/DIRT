import time
import traceback

# import mahotas as m
import numpy as np
from graph_tool import Graph
import graph_tool.topology as gt
from graph_tool.libgraph_tool_core import VertexBase
from typing import List

from utils import compare_lists


def extract_graph(image: np.ndarray, diameter, x_scale: float, y_scale: float) -> Graph:
    print('Building graph'),
    start = time.time()
    graph = Graph(directed=False)
    added_vertices_line2 = []
    sum_add_vertices = 0
    v_list_line2 = []
    percent_old = 0
    counter = 0
    vertex_prop = graph.new_vertex_property('object')
    edge_prop = graph.new_edge_property('object')
    edge_prop_w = graph.new_edge_property("float")
    height, width = np.shape(image)

    if x_scale > 0 and y_scale > 0:
        avg_scale = (x_scale + y_scale) / 2
    else:
        avg_scale = 1.
        x_scale = 1.
        y_scale = 1.

    # sweep over each line in the image except the last line
    for idx, i in enumerate(image[:len(image) - 2]):
        # get foreground indices in the current line of the image and make vertices
        counter += 1
        percent = (float(counter) / float(height)) * 100
        if percent_old + 10 < percent:
            print(str(np.round(percent, 1)) + '% '),
            percent_old = percent

        line1 = np.where(i)
        if len(line1[0]) > 0:
            line1 = set(line1[0]).difference(set(added_vertices_line2))
            vl = graph.add_vertex(len(list(line1)))

            if len(line1) > 1:
                v_list = v_list_line2 + list(vl)
            else:
                v_list = v_list_line2 + [vl]
            line1 = added_vertices_line2 + list(line1)
            for jdx, j in enumerate(line1):
                vertex_prop[v_list[jdx]] = {'imgIdx': (j, idx), 'coord': (float(j) * x_scale, float(idx) * y_scale),
                                            'abs_coord': (float(j), float(idx)),
                                            'nrOfPaths': 0, 'diameter': float(diameter[idx][j]) * avg_scale}

            # keep order of the inserted vertices
            sum_add_vertices += len(line1)
            added_vertices_line2 = []
            v_list_line2 = []

            # connect foreground indices to neighbours in the next line
            for v1 in line1:
                va = v_list[line1.index(v1)]
                diagonal_left = diagonal_right = True
                try:
                    if image[idx][v1 - 1]:
                        diagonal_left = False
                        vb = v_list[line1.index(v1 - 1)]
                        e = graph.add_edge(va, vb)
                        edge_prop[e] = {'coord1': vertex_prop[va]['coord'], 'coord2': vertex_prop[vb]['coord'],
                                        'w': ((vertex_prop[va]['diameter'] + vertex_prop[vb]['diameter']) / 2),
                                        'RTP': False}
                        edge_prop_w[e] = 2. / (edge_prop[e]['w'] ** 2)
                except:
                    print('Boundary vertex at: ' + str([v1, idx - 1]) + ' image size: ' + str([width, height]))
                    pass

                try:
                    if image[idx][v1 + 1]:
                        diagonal_right = False
                        vb = v_list[line1.index(v1 + 1)]
                        e = graph.add_edge(va, vb)
                        edge_prop[e] = {'coord1': vertex_prop[va]['coord'], 'coord2': vertex_prop[vb]['coord'],
                                        'w': ((vertex_prop[va]['diameter'] + vertex_prop[vb]['diameter']) / 2),
                                        'RTP': False}
                        edge_prop_w[e] = 2. / (edge_prop[e]['w'] ** 2)
                except:
                    print('Boundary vertex at: ' + str([v1 + 1, idx]) + ' image size: ' + str([width, height]))
                    pass  # just if we are out of bounds

                try:
                    if image[idx + 1][v1]:
                        diagonal_right = False
                        diagonal_left = False
                        vNew = graph.add_vertex()
                        vertex_prop[vNew] = {'imgIdx': (v1, idx + 1),
                                             'coord': (float(v1) * x_scale, float(idx + 1) * y_scale), 'nrOfPaths': 0,
                                             'abs_coord': (float(j), float(idx)),
                                             'diameter': float(diameter[idx + 1][v1]) * avg_scale}
                        v_list_line2.append(vNew)
                        e = graph.add_edge(v_list[line1.index(v1)], vNew)
                        edge_prop[e] = {'coord1': vertex_prop[va]['coord'], 'coord2': vertex_prop[vNew]['coord'],
                                        'w': ((vertex_prop[va]['diameter'] + vertex_prop[vNew]['diameter']) / 2),
                                        'RTP': False}
                        edge_prop_w[e] = 1. / (edge_prop[e]['w'] ** 2)
                        if v1 not in added_vertices_line2: added_vertices_line2.append(v1)
                except:
                    print('Boundary vertex at: ' + str([v1, idx + 1]) + ' image size: ' + str([width, height]))
                    pass

                try:
                    if diagonal_right and image[idx + 1][v1 + 1]:
                        vNew = graph.add_vertex()
                        vertex_prop[vNew] = {'imgIdx': (v1 + 1, idx + 1),
                                             'coord': (float(v1 + 1) * x_scale, float(idx + 1) * y_scale),
                                             'abs_coord': (float(j), float(idx)),
                                             'nrOfPaths': 0,
                                             'diameter': float(diameter[idx + 1][v1 + 1]) * avg_scale}
                        v_list_line2.append(vNew)
                        e = graph.add_edge(v_list[line1.index(v1)], vNew)
                        edge_prop[e] = {'coord1': vertex_prop[va]['coord'], 'coord2': vertex_prop[vNew]['coord'],
                                        'w': ((vertex_prop[va]['diameter'] + vertex_prop[vNew]['diameter']) / 2),
                                        'RTP': False}
                        edge_prop_w[e] = 1.41 / (edge_prop[e]['w'] ** 2)
                        if v1 + 1 not in added_vertices_line2: added_vertices_line2.append(v1 + 1)
                except:
                    print('Boundary vertex at: ' + str([v1 + 1, idx + 1]) + ' image size: ' + str([width, height]))
                    pass

                try:
                    if diagonal_left and image[idx + 1][v1 - 1]:
                        vNew = graph.add_vertex()
                        vertex_prop[vNew] = {'imgIdx': (v1 - 1, idx + 1),
                                             'coord': (float(v1 - 1) * x_scale, float(idx + 1) * y_scale),
                                             'abs_coord': (float(j), float(idx)),
                                             'nrOfPaths': 0,
                                             'diameter': float(diameter[idx + 1][v1 - 1]) * avg_scale}
                        v_list_line2.append(vNew)
                        e = graph.add_edge(v_list[line1.index(v1)], vNew)
                        edge_prop[e] = {'coord1': vertex_prop[va]['coord'], 'coord2': vertex_prop[vNew]['coord'],
                                        'w': ((vertex_prop[va]['diameter'] + vertex_prop[vNew]['diameter']) / 2),
                                        'RTP': False}
                        edge_prop_w[e] = 1.41 / (edge_prop[e]['w'] ** 2)
                        if v1 - 1 not in added_vertices_line2: added_vertices_line2.append(v1 - 1)
                except:
                    print('Boundary vertex at: ' + str([v1 - 1, idx + 1]) + ' image size: ' + str([width, height]))
                    pass

                try:
                    if not image[idx][v1 + 1] and \
                            not image[idx][v1 - 1] and \
                            not image[idx + 1][v1] and \
                            not diagonal_left and \
                            not diagonal_right:
                        print('tip detected')
                        if not image[idx - 1][v1 - 1] and \
                                not image[idx - 1][v1 + 1] and \
                                not image[idx - 1][v1]:
                            print('floating pixel')
                except:
                    pass

    graph.edge_properties["ep"] = edge_prop
    graph.edge_properties["w"] = edge_prop_w
    graph.vertex_properties["vp"] = vertex_prop
    largest_comp = gt.GraphView(graph, vfilt=gt.label_largest_component(graph))
    print(f"Complete graph with {graph.num_vertices()} vertices built in {time.time() - start}")

    if largest_comp.num_vertices() != graph.num_vertices():
        largest_comp_reduction = round(
            float((graph.num_vertices() - largest_comp.num_vertices())) / float(graph.num_vertices()), 2)
        print(
            f"Largest component has {largest_comp.num_vertices()} vertices ({largest_comp_reduction}% vertices omitted)")

    return largest_comp


def initial_branching_point(graph: Graph, image_height: int) -> VertexBase:
    print(f"Finding vertex at y = {image_height}")
    height = image_height
    target = None
    vp = graph.vertex_properties["vp"]

    for vertex in graph.vertices():
        count = 0
        for _ in vertex.out_neighbours():
            count += 1
            if count > 2:
                break
        if count > 2:
            if vp[vertex]['imgIdx'][1] < height:
                height = vp[vertex]['imgIdx'][1]
                target = vertex

    return target


def lateral(graph: Graph, image_height: int) -> VertexBase:
    print('Finding root vertex lateral X')
    target = None
    height = image_height
    vp = graph.vertex_properties["vp"]

    for vertex in graph.vertices():
        if vp[vertex]['imgIdx'][1] < height:
            height = vp[vertex]['imgIdx'][1]
            target = vertex

    return target


def last_branch_point(graph: Graph) -> VertexBase:
    path = 0
    target = None
    vp = graph.vertex_properties["vp"]

    for vertex in graph.vertices():
        try:
            if vp[vertex]['imgIdx'][1] > path:
                path = vp[vertex]['imgIdx'][1]
                target = vertex
        except:
            pass

    return target


def thickest_full_path(graph: Graph, branch_point: VertexBase) -> List[VertexBase]:
    print('Finding thickest diameter trace path')
    weights = graph.edge_properties["w"]
    path = []

    if graph.num_vertices() > 0:
        path, _ = gt.shortest_path(
            graph,
            branch_point,
            last_branch_point(graph),
            weights=weights,
            pred_map=None)

    return path


def hypocotol_cluster(graph: Graph, path: List[VertexBase]) -> (List[VertexBase], List[float]):
    print('Find hypocotol cluster')
    radius = []
    branching_paths = []
    branching_vertices = []
    vp = graph.vertex_properties["vp"]

    for vertex in path:
        branching_paths.append(vp[vertex]['nrOfPaths'])
        branching_vertices.append(vertex)

    for vertex in branching_vertices:
        radius.append(vp[vertex]['diameter'])

    bp = []
    rad = []
    temp_avg = 0.
    counter = 0.
    for vertex in range(len(branching_vertices) - 1):
        if branching_paths[vertex] == branching_paths[vertex + 1]:
            temp_avg += radius[vertex]
            counter += 1
        elif counter > 0:
            temp_avg = temp_avg / counter
            rad.append(temp_avg)
            bp.append(branching_paths[vertex])
            counter = 0.
            temp_avg = 0.

    return bp, rad


def cluster_overlay_image(
        graph: Graph,
        path: List[VertexBase],
        image: np.ndarray,
        x_scale: float,
        y_scale: float,
        c1x,
        c1y,
        c2x,
        c2y,
        c3x=None,
        c3y=None):
    print('Creating cluster overlay image')
    # image = m.as_rgb(image, image, image)
    vp = graph.vertex_properties["vp"]

    for vertex in path:
        if vp[vertex]['nrOfPaths'] in c1y:
            y = int(vp[vertex]['imgIdx'][0])
            x = int(vp[vertex]['imgIdx'][1])
            try:
                image[x][y] = (125, 0, 0)
            except:
                pass

            diameter = vp[vertex]['diameter'] / (x_scale / 2 + y_scale / 2)
            diameter = diameter * 1.5

            for i in range(int(diameter)):
                try:
                    image[x][y + i] = (125, 0, 0)
                except:
                    pass
                try:
                    image[x][y - i] = (125, 0, 0)
                except:
                    pass
                try:
                    image[x - i][y] = (125, 0, 0)
                except:
                    pass
                try:
                    image[x + i][y] = (125, 0, 0)
                except:
                    pass
        elif vp[vertex]['nrOfPaths'] in c2y:
            y = int(vp[vertex]['imgIdx'][0])
            x = int(vp[vertex]['imgIdx'][1])
            try:
                image[x][y] = (125, 0, 0)
            except:
                pass

            diameter = vp[vertex]['diameter'] / (x_scale / 2 + y_scale / 2)
            diameter = diameter * 1.5

            for i in range(int(diameter)):
                try:
                    image[x][y + i] = (0, 125, 0)
                except:
                    pass
                try:
                    image[x][y - i] = (0, 125, 0)
                except:
                    pass
                try:
                    image[x - i][y] = (0, 125, 0)
                except:
                    pass
                try:
                    image[x + i][y] = (0, 125, 0)
                except:
                    pass
            y = int(vp[vertex]['imgIdx'][0])
            x = int(vp[vertex]['imgIdx'][1])
            try:
                image[x][y] = (0, 0, 125)
            except:
                pass

            diameter = vp[vertex]['diameter'] / (x_scale / 2 + y_scale / 2)
            diameter = diameter * 1.5

            for i in range(int(diameter)):
                try:
                    image[x][y + i] = (0, 0, 125)
                except:
                    pass
                try:
                    image[x][y - i] = (0, 0, 125)
                except:
                    pass
                try:
                    image[x - i][y] = (0, 0, 125)
                except:
                    pass
                try:
                    image[x + i][y] = (0, 0, 125)
                except:
                    pass

    return image


def lateral_lengths(graph: Graph, tips: List[VertexBase], scale: float):
    branch_points = []
    lateral_lengths = []
    vp = graph.vertex_properties["vp"]

    for tip in tips:
        neighbors = list(tip.all_neighbors())
        if len(neighbors) > 0:
            for neighbor in neighbors:
                d = float(vp[graph.vertex(neighbor)]['diameter'])
                radius = int(d / (1.0 if scale == 0.0 else scale))  # convert radius at branching point to pixels
                if radius > 0:
                    break

        # remove the radius of the main trunk from the lateral length
        # to obtain the emerging lateral length from the surface
        if radius + 2 < len(neighbors):
            lateral_lengths.append(tip[radius:])
            branch_points.append(tip[0])

    return lateral_lengths, branch_points


def root_tips(graph: Graph, path: List[VertexBase]):
    print('Finding root tips')
    tips = []
    diameters = []
    heights = []
    widths = []
    vp = graph.vertex_properties["vp"]

    for vertex in graph.vertices():
        count = 0
        for _ in vertex.out_neighbours():
            count += 1
        if count == 3:
            widths.append(vp[vertex]['coord'][1])
        if count <= 1:
            if vertex != path[len(path) - 1]:
                tips.append(vertex)
                diameters.append(vp[vertex]['diameter'])
                heights.append(vp[vertex]['coord'][0])

    percent_90 = np.max(heights) * 0.9
    tmp_dia_90 = []
    for vertex in list(np.where(heights >= percent_90)[0]):
        tmp_dia_90.append(diameters[vertex])

    # dia90 = np.max(tmp_dia_90) if tmp_dia_90 else 0   TODO why is this line here?
    dia90 = np.max(diameters) if diameters else 0
    max_depth = np.max(heights) if heights else 0
    max_width = np.max(widths) - np.min(widths) if widths else 0

    return tips, diameters, dia90, max_depth, max_width, heights


def root_tip_paths(graph: Graph, path: List[VertexBase]):
    print('Finding root tip paths')
    tips, dia, dia90, root_depth, root_width, tip_height = root_tips(graph, path)
    paths = []

    if tips == -1:
        print('ERROR: No root tips found')
        return -1
    try:
        tips.remove(path[0])
    except:
        pass

    vp = graph.vertex_properties["vp"]
    ep = graph.edge_properties["ep"]
    weights = graph.edge_properties["w"]
    percent_old = 0
    for idx, vertex in enumerate(tips):
        percent = (float(idx) / float(len(tips))) * 100
        if percent_old + 5 < percent:
            print(str(np.round(percent, 1)) + '% '),
            percent_old = percent
        try:
            path, edges = gt.shortest_path(
                graph,
                path[0],
                vertex,
                weights=weights,
                pred_map=None)
            paths_temp = []
            for e in path:
                paths_temp.append(e)
            split = compare_lists(path, paths_temp)
            paths.append(paths_temp[split:])

            for e in reversed(path):
                vp[e]['nrOfPaths'] += 1
            for e in edges:
                ep[e]['RTP'] = True
        except:
            print(f"No Dijkstra path for tip {idx}: {traceback.format_exc()}")
            pass
    print(f"Found {len(paths)} root tip paths")

    return paths, tips, dia, dia90, root_depth, root_width, tip_height


def root_tip_skeleton(graph: Graph, path: List[VertexBase]):
    paths, tips, diameter, diameter_90, depth, width, height = root_tip_paths(graph, path)
    graph_copy = graph.copy()
    ep = graph.edge_properties["ep"]

    for edge in graph.edges():
        if not ep[edge]['RTP']:
            graph_copy.remove_edge(edge)

    return graph_copy, paths, tips, diameter, diameter_90, depth, width, height
