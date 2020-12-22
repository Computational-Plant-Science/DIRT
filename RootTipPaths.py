'''
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

'''
import sys
import traceback

'''
# external library imports
'''
import numpy as np
import scipy.optimize as sp
import graph_tool.topology as gt

'''
# standard python imports
'''
import time


class RootTipPaths(object):
    '''
    classdocs
    '''

    def __init__(self, io):
        '''
        Constructor
        '''
        self.__RTP = []
        self.__io = io
        self.__id = io.getID()
        self.__currentIdx = io.getCurrentID()
        self.__medianTipDiameter = 0.0
        self.__meanTipDiameter = 0.0
        self.__90TipDiameter = 0.0
        self.__expF = 0.0
        self.__tips = []
        self.__rootingDepth = 0.0
        self.__rootWidth = 0.0

    def compareTwoOrderedLists(self, l1, l2):

        if len(l1) > len(l2): l1 = l1[:len(l2)]
        if len(l1) < len(l2): l2 = l2[:len(l1)]
        half = len(l1) // 2
        split = 0

        if l1[half] == l2[half]:
            if half == 0:
                return 0
            elif half == len(l1) - 1:
                return len(l1) - 1
            elif l1[half + 1] == l2[half + 1]:
                split = self.compareTwoOrderedLists(l1[half:], l2[half:])
            else:
                split = 0
            return half + split
        else:
            if half == 0: return -1
            split = self.compareTwoOrderedLists(l1[:half], l2[:half])
            return split

    def getAllTips(self):
        return self.__tips

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

    def getRootTipPaths(self, thickestPath, G):
        print('Calculating Root-Tip Paths')

        CPVIDX = []
        for i in thickestPath:
            CPVIDX.append(G.vertex_index[i])
        vprop = G.vertex_properties["vp"]
        eprop = G.edge_properties["ep"]
        epropW = G.edge_properties["w"]
        if len(self.__RTP) == 0:
            tips = self.getTips(thickestPath, G)
            RTP = []
            # print '***** TIPS VAR ******'
            # print tips
            if tips == -1:
                print('ERROR: No tips found')
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
                    split = self.compareTwoOrderedLists(CPVIDX, RTPTmp)
                    RTP.append(RTPTmp[split:])

                    for j in reversed(path):
                        vprop[j]['nrOfPaths'] += 1
                    for j in edges:
                        eprop[j]['RTP'] = True
                except:
                    print(traceback.format_exc())
                    print('ERROR: in def getRootTipPaths(self,thickestPath,G): no dijkstra path at ' + str(
                        idx) + ' in tips')
                    pass
            print('Number of Root-Tip Paths: ' + str(len(RTP)))
            self.__RTP = RTP
        print('RTP done!')
        return RTP, tips

    def getTips(self, thickestPath, G, counter=None):
        tips = []
        tipDia = []
        tipHeight = []
        rootW = []
        vprop = G.vertex_properties["vp"]
        for i in G.vertices():
            count = 0
            for _ in i.out_neighbours():
                count += 1
            if count == 3:
                rootW.append(vprop[i]['coord'][1])
            if count <= 1:
                if i != thickestPath[len(thickestPath) - 1]:
                    tips.append(G.vertex_index[i])
                    tipDia.append(vprop[i]['diameter'])
                    tipHeight.append(vprop[i]['coord'][0])
        self.__medianTipDiameter = np.median(tipDia)
        print('Median Tip Diameter: ' + str(self.__medianTipDiameter))
        self.__meanTipDiameter = np.mean(tipDia)
        print('Mean Tip Diameter: ' + str(self.__meanTipDiameter))
        self.__io.saveArray(tipDia, self.__io.getHomePath() + 'Plots/' + self.__io.getFileName() + '_TipDiaHisto')
        self.__io.saveArray(tipDia, self.__io.getHomePath() + 'Plots/' + self.__io.getFileName() + '_TipDiaHeightX')
        self.__io.saveArray(tipHeight, self.__io.getHomePath() + 'Plots/' + self.__io.getFileName() + '_TipDiaHeightY')

        try:
            percent90 = np.max(tipHeight) * 0.9
            idx = list(np.where(tipHeight >= percent90)[0])
            tmpdia90 = []
            for i in idx:
                tmpdia90.append(tipDia[i])

            if tmpdia90:
                dia90 = np.max(tmpdia90)
            else:
                dia90 = 0
            self.__90TipDiameter = dia90
            if tipDia:
                self.__90TipDiameter = np.max(tipDia)
            else:
                self.__90TipDiameter = 0
            if tipHeight:
                self.__rootingDepth = np.max(tipHeight)
            else:
                self.__rootingDepth = 0
            if rootW:
                self.__rootWidth = np.max(rootW) - np.min(rootW)
            else:
                self.__rootWidth = 0
        except:
            pass
        return tips

    def getRTPSkeleton(self, thickestPath, G, newRTp=False):
        eprop = G.edge_properties["ep"]
        if newRTp == True: self.__RTP = []
        if len(self.__RTP) == 0:
            startT = time.time()
            RTP, tips = self.getRootTipPaths(thickestPath, G)
            self.__RTP = RTP
            print('RTPs computed in ' + str(time.time() - startT) + 's')
        print('calculating RTP Skeleton')

        rtpSkel = G.copy()
        for e in G.edges():
            if eprop[e]['RTP'] == False:
                rtpSkel.remove_edge(e)

        return rtpSkel, len(
            self.__RTP), self.__medianTipDiameter, self.__meanTipDiameter, self.__90TipDiameter, self.__RTP, tips, self.__rootingDepth, self.__rootWidth
