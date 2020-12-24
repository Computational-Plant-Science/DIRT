"""
RootTipPaths.py

This module is used to compute the RTPs and the RTP skeleton as described in the paper.

The code is free for non-commercial use.
Please contact the author for commercial use.

Please cite the DIRT Paper if you use the code for your scientific project.

Bucksch et al., 2014 "Image-based high-throughput field phenotyping of crop roots", Plant Physiology

-------------------------------------------------------------------------------------------
Author: Alexander Bucksch
School of Biology and Interactive computing
Georgia Institute of Technology

Mail: bucksch@gatech.edu
Web: http://www.bucksch.nl
-------------------------------------------------------------------------------------------

Copyright (c) 2014 Alexander Bucksch
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided
    with the distribution.

  * Neither the name of the DIRT Developers nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import time
import traceback

import graph_tool.topology as gt
import numpy as np
import scipy.optimize as sp

from options import DIRTOptions


class RootTips(object):
    def __init__(self, options: DIRTOptions):
        self.__RTP = []
        self.__medianTipDiameter = 0.0
        self.__meanTipDiameter = 0.0
        self.__90TipDiameter = 0.0
        self.__expF = 0.0
        self.__tips = []
        self.__rootingDepth = 0.0
        self.__rootWidth = 0.0
        self.__options = options

    def __compare_lists(self, l1, l2):

        if len(l1) > len(l2): l1 = l1[:len(l2)]
        if len(l1) < len(l2): l2 = l2[:len(l1)]
        half = len(l1) // 2

        if l1[half] == l2[half]:
            if half == 0:
                return 0
            elif half == len(l1) - 1:
                return len(l1) - 1
            elif l1[half + 1] == l2[half + 1]:
                split = self.__compare_lists(l1[half:], l2[half:])
            else:
                split = 0
            return half + split
        else:
            if half == 0: return -1
            split = self.__compare_lists(l1[:half], l2[:half])
            return split

    def model_func(self, t, A, K, C):
        return A * np.exp(K * t) + C

    def model_func_dia(self, t, A, K, C, dia):
        return A * dia ** (K * t) + C

    def fit_exp_nonlinear(self, t, y, dia):
        opt_parms, _ = sp.curve_fit(self.model_func, t, y, maxfev=100000)
        A, K, C = opt_parms
        fit_y = self.model_func(t, A, K, C)
        return fit_y, A, K, C

    def fit_exp_linear(self, t, y, C, dia):
        y = np.array(y) - C
        y = np.log(y)
        K, A_log = np.polyfit(t, y, 1)
        A = np.exp(A_log)
        fit_y = self.model_func(np.array(t), A, K, C)
        return fit_y, A, K, C

    def get_paths(self, thickestPath, G):
        print('Finding root tip paths')
        CPVIDX = []
        for i in thickestPath:
            CPVIDX.append(G.vertex_index[i])
        vprop = G.vertex_properties["vp"]
        eprop = G.edge_properties["ep"]
        epropW = G.edge_properties["w"]
        if len(self.__RTP) == 0:
            tips = self.get_tips(thickestPath, G)
            RTP = []
            if tips == -1:
                print('ERROR: No root tips found')
                return -1
            try:
                tips.remove(G.vertex_index(thickestPath[0]))
            except:
                pass

            percentOld = 0
            for idx, i in enumerate(tips):

                percent = (float(idx) / float(len(tips))) * 100
                if percentOld + 5 < percent:
                    print(str(np.round(percent, 1)) + '% '),
                    percentOld = percent
                try:
                    path, edges = gt.shortest_path(G, thickestPath[0], G.vertex(i), weights=epropW, pred_map=None)
                    RTPTmp = []
                    for k in path:
                        RTPTmp.append(G.vertex_index[k])
                    split = self.__compare_lists(CPVIDX, RTPTmp)
                    RTP.append(RTPTmp[split:])

                    for j in reversed(path):
                        vprop[j]['nrOfPaths'] += 1
                    for j in edges:
                        eprop[j]['RTP'] = True
                except:
                    print(f"No Dijkstra path for tip {idx}: {traceback.format_exc()}")
                    pass
            print(f"Found {len(RTP)} root-tip paths")
            self.__RTP = RTP
        return RTP, tips

    def get_tips(self, thickest_path, graph):
        tips = []
        tip_dia = []
        tip_height = []
        root_w = []
        vprop = graph.vertex_properties["vp"]
        for i in graph.vertices():
            count = 0
            for _ in i.out_neighbours():
                count += 1
            if count == 3:
                root_w.append(vprop[i]['coord'][1])
            if count <= 1:
                if i != thickest_path[len(thickest_path) - 1]:
                    tips.append(graph.vertex_index[i])
                    tip_dia.append(vprop[i]['diameter'])
                    tip_height.append(vprop[i]['coord'][0])
        self.__medianTipDiameter = np.median(tip_dia)
        print('Median tip diameter: ' + str(self.__medianTipDiameter))
        self.__meanTipDiameter = np.mean(tip_dia)
        print('Mean tip diameter: ' + str(self.__meanTipDiameter))
        np.savetxt(f"{self.__options.input_stem}.tip.dia.histo.txt", tip_dia, delimiter=',')
        np.savetxt(f"{self.__options.input_stem}.tip.dia.height.x.txt", tip_dia, delimiter=',')
        np.savetxt(f"{self.__options.input_stem}.tip.dia.height.y.txt", tip_height, delimiter=',')

        try:
            percent_90 = np.max(tip_height) * 0.9
            idx = list(np.where(tip_height >= percent_90)[0])
            tmp_dia_90 = []
            for i in idx:
                tmp_dia_90.append(tip_dia[i])

            if tmp_dia_90:
                dia90 = np.max(tmp_dia_90)
            else:
                dia90 = 0
            self.__90TipDiameter = dia90
            if tip_dia:
                self.__90TipDiameter = np.max(tip_dia)
            else:
                self.__90TipDiameter = 0
            if tip_height:
                self.__rootingDepth = np.max(tip_height)
            else:
                self.__rootingDepth = 0
            if root_w:
                self.__rootWidth = np.max(root_w) - np.min(root_w)
            else:
                self.__rootWidth = 0
        except:
            print(traceback.format_exc())
        return tips

    def get_rtp_skeleton(self, thickest_path, graph, new_rtp=False):
        eprop = graph.edge_properties["ep"]
        if new_rtp: self.__RTP = []
        if len(self.__RTP) == 0:
            start = time.time()
            paths, tips = self.get_paths(thickest_path, graph)
            self.__RTP = paths
            print('RTPs computed in ' + str(time.time() - start) + 's')
        print('Finding root tip path skeleton')

        rtpSkel = graph.copy()
        for e in graph.edges():
            if not eprop[e]['RTP']:
                rtpSkel.remove_edge(e)

        return rtpSkel, len(
            self.__RTP), self.__medianTipDiameter, self.__meanTipDiameter, self.__90TipDiameter, self.__RTP, tips, self.__rootingDepth, self.__rootWidth
