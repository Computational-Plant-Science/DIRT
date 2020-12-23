#! /nv/hp10/adas30/bin/python
'''
----------------------------------------------------------------------------------------------------
DIRT 1.1 - An automatic high throughput root phenotyping platform
Web interface by Abhiram Das - adas30@biology.gatech.edu

http://dirt.iplantcollaborative.org

University of Georgia

The software is written in:
- python 3 (https://www.python.org)

The software depends on:
- the graphtools package (http://graph-tool.skewed.de) 
- the mahotas package (http://luispedro.org/software/mahotas)
- the numpy package (http://sourceforge.net/projects/numpy/)
- the scipy package (http://www.scipy.org/SciPy)

Optionally binaries of can be used for tag recognition:

- tesseract (https://code.google.com/p/tesseract-ocr/)
- zbar (http://zbar.sourceforge.net)

The software uses free code that had no references when found on the net:
- http://www.daniweb.com/software-development/python/threads/31449/k-means-clustering

The software uses modified code from the scikit.image:
- adaptive thresholding in Masking.py (http://scikit-image.org)

The software uses modified code from Kyle Fox:
- fixOrientation.py: https://github.com/kylefox/python-image-orientation-patch


Please cite the DIRT Paper if you use the code for your scientific project.

Bucksch et al., 2014 "Image-based high-throughput field phenotyping of crop roots", Plant Physiology


----------------------------------------------------------------------------------------------------
Author: Alexander Bucksch
Department of Plant Biology
Warnell School of Forestry and Natural Resources
Institute of Bioinformatics
University of Georgia

Mail: bucksch@uga.edu
Web: http://www.computational-plant-science.org
----------------------------------------------------------------------------------------------------

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

'''
# external library imports
'''
import scipy
import imageio
import skimage

'''
# internal library imports
'''
import dirtIO
import Segmentation
import Preprocessing
import Skeleton
import Analysis
import RootTipPaths
from fixImageOrientation import *

'''
# standard python imports
'''
import glob
import os
import pickle
import csv
import sys
import time
from collections import OrderedDict

'''
#global defs
'''
allCrown = []
allPara = []
f = []
imgID = None
io = dirtIO.IO()
options = []
scale = None
ID = None
stemCorrection = False
maxExRoot = None
traitDict = OrderedDict()


def init(fpath, io):
    oldpath = os.getcwd()
    # io.setHomePath(fpath)
    # if not os.path.exists(fpath):
    #     os.mkdir(fpath)
    # os.chdir(fpath)
    # print(os.getcwd())
    io.setServerPath(os.getcwd())
    # if not os.path.exists('tmp'):
    #     os.mkdir('tmp')

    os.chdir(oldpath)
    readTraits(options[12][1])


def readTraits(myFilePath='./traits.csv'):
    global traitDict
    print("TRAITS DIRECTORY: " + os.getcwd())
    # check to make sure its a file not a sub folder
    if (os.path.isfile(myFilePath) and myFilePath.endswith(".csv")):
        with open(myFilePath, 'U') as csvfile:
            # sniff to find the format
            fileDialect = csv.Sniffer().sniff(csvfile.read(1024))
            csvfile.seek(0)
            # read the CSV file into a dictionary
            dictReader = csv.reader(csvfile, dialect=fileDialect)
            for row in dictReader:
                try:
                    traitDict[row[0]] = bool(int(row[1]))
                except:
                    print('invalid entry in trait file: ' + str(row))
                    pass
    print(traitDict)
    return


def readOptions():
    global options
    if len(sys.argv) == 13:
        options.append([0, os.path.dirname(sys.argv[1])])
        options.append([0, os.path.basename(sys.argv[1])])
        options.append([0, sys.argv[2]])
        options.append([0, sys.argv[3]])
        options.append([0, sys.argv[4]])
        options.append([0, sys.argv[5]])
        options.append([0, sys.argv[6]])
        options.append([0, sys.argv[7]])
        options.append([0, sys.argv[8]])
        options.append([0, sys.argv[9]])
        options.append([0, sys.argv[10]])
        options.append([0, sys.argv[11]])
        options.append([0, sys.argv[12]])

    else:
        with open('./options.csv', 'U') as csvfile:
            filedata = csv.reader(csvfile)
            for i in filedata:
                options.append(i)

    return options


def ifAnyKeyIsTrue(listOfKeys):
    for i in listOfKeys:
        if traitDict[i] == True:
            return True
    return False


def threadSegmentation(filepath, imgFile, imgID, maxExRoot, rootCrown, marker):
    global io
    global scale
    global stemCorrection

    stemCorrection = bool(int(options[8][1]))
    io.setFileName(imgFile)
    io.setidIdx(imgID)
    prep = Preprocessing.Preprocessing(io)
    print('segmenting file: ' + imgFile + '\n')

    image_file_path = os.path.join(options[0][1], imgFile)

    if os.path.isfile(image_file_path):
        # fix orientation of the image in tiff and Jpg files
        fix_orientation(image_file_path, save_over=True)
        img = imageio.imread(image_file_path, as_gray=True)

    else:
        print('Image not readable')
        img = []

    if len(img) > 0:
        currT = time.time()
        Failed, tagExtract, circleRatio, circleWidth, circleHeight = prep.prepocess(img, rootCrown,
                                                                                    scale=float(options[3][1]),
                                                                                    nrExRoot=maxExRoot, marker=marker,
                                                                                    stemCorrection=stemCorrection)
        print('Segmentation finished in ' + str(time.time() - currT) + 's')
        if Failed == False:
            xScale = scale / float(circleWidth)
            yScale = scale / float(circleHeight)
            if xScale <= 0.0: xScale = 1.
            if yScale <= 0.0: yScale = 1.
            para = [int(imgID), io.getFileName(), Failed, tagExtract, circleRatio, circleWidth, circleHeight, xScale,
                    yScale, -1, -1]
            if maxExRoot > 1:
                for _ in range(maxExRoot):
                    allPara.append(para)
            else:
                allPara.append(para)
        else:
            xScale = scale / float(1.0)
            yScale = scale / float(1.0)
            circleRatio = 1.0
            circleWidth = 1.0
            circleHeight = 1.0
            para = [int(imgID), io.getFileName(), Failed, tagExtract, circleRatio, circleWidth, circleHeight, xScale,
                    yScale, -1, -1]
            if maxExRoot > 1:
                for _ in range(maxExRoot):
                    allPara.append(para)
            else:
                allPara.append(para)


def threadCrown(filepath):
    global io

    rtpSkel = -1
    crownT = OrderedDict()
    imgL = []
    stemCorrection = bool(int(options[8][1]))

    print(io.getHomePath())
    oldHome = io.getHomePath()
    os.chdir(io.getHomePath())
    f = glob.glob('*.crown.png')
    for (counter, i) in enumerate(f):
        io.setFileName(os.path.basename(i))
        io.setidIdx(imgID)

        print('processing Crown file: ' + i)
        xScale = allPara[counter][7]
        yScale = allPara[counter][8]
        analysis = Analysis.Analysis(io, (xScale + yScale) / 2)
        rtp = RootTipPaths.RootTipPaths(io)

        try:
            img = imageio.imread(i, as_gray=True)
        except:
            print('Image not readable')
            img = -1

        if len(img) > 0:
            seg = Segmentation.Segmentation(img, io)
            imgL = seg.label()
            print('compute root profile')
            currT = time.time()
            if ifAnyKeyIsTrue(
                    ['AVG_DENSITY', 'WIDTH_MED', 'WIDTH_MAX', 'DIA_STM_SIMPLE', 'D10', 'D20', 'D30', 'D40', 'D50',
                     'D60', 'D70', 'D80', 'D90', 'DS10', 'DS20', 'DS30', 'DS40', 'DS50', 'DS60', 'DS70', 'DS80', 'DS90',
                     'AREA', 'ANG_TOP', 'ANG_BTM']):
                crownT['AVG_DENSITY'], crownT['WIDTH_MED'], crownT['WIDTH_MAX'], crownT['D10'], crownT['D20'], crownT[
                    'D30'], crownT['D40'], crownT['D50'], crownT['D60'], crownT['D70'], crownT['D80'], crownT['D90'], \
                crownT['DS10'], crownT['DS20'], crownT['DS30'], crownT['DS40'], crownT['DS50'], crownT['DS60'], crownT[
                    'DS70'], crownT['DS80'], crownT['DS90'], crownT['AREA'], crownT['DIA_STM_SIMPLE'], crownT[
                    'ANG_TOP'], crownT['ANG_BTM'] = analysis.getWidthOverHeight(imgL, xScale, yScale)
                print('Mask traits computed ' + str(time.time() - currT) + 's')

            if ifAnyKeyIsTrue(
                    ['DIA_STM', 'TD_MED', 'TD_AVG', 'STA_RANGE', 'STA_DOM_I', 'STA_DOM_II', 'STA_25_I', 'STA_25_II',
                     'STA_50_I', 'STA_50_II', 'STA_75_I', 'STA_75_II', 'STA_90_I', 'STA_90_II', 'RTA_DOM_I',
                     'RTA_DOM_II', 'STA_MIN', 'STA_MAX', 'STA_MED', 'RTA_RANGE', 'RTA_MIN', 'RTA_MAX', 'RTA_MED',
                     'NR_RTP_SEG_I', 'NR_RTP_SEG_II', 'ADVT_COUNT', 'BASAL_COUNT', 'ADVT_ANG', 'BASAL_ANG', 'HYP_DIA',
                     'TAP_DIA', 'MAX_DIA_90', 'DROP_50', 'CP_DIA25', 'CP_DIA50', 'CP_DIA75', 'CP_DIA90', 'SKL_DEPTH',
                     'SKL_WIDTH']):
                currT = time.time()
                skel = Skeleton.Skeleton(imgL)
                testSkel, testDia = skel.skel(imgL)
                imageio.imwrite(io.getFileName() + '.skel.png', skimage.img_as_uint(testSkel))
                print('Medial axis computed ' + str(time.time() - currT) + 's')
                currT = time.time()
                path, skelGraph, crownT['DIA_STM'], skelSize = seg.findThickestPath(testSkel, testDia, xScale, yScale)
                allPara[counter][10] = skelSize
                print('Central path computed ' + str(time.time() - currT) + 's')

            if ifAnyKeyIsTrue(
                    ['TD_MED', 'TD_AVG', 'STA_RANGE', 'STA_DOM_I', 'STA_DOM_II', 'STA_25_I', 'STA_25_II', 'STA_50_I',
                     'STA_50_II', 'STA_75_I', 'STA_75_II', 'STA_90_I', 'STA_90_II', 'RTA_DOM_I', 'RTA_DOM_II',
                     'STA_MIN', 'STA_MAX', 'STA_MED', 'RTA_RANGE', 'RTA_MIN', 'RTA_MAX', 'RTA_MED', 'NR_RTP_SEG_I',
                     'NR_RTP_SEG_II', 'ADVT_COUNT', 'BASAL_COUNT', 'ADVT_ANG', 'BASAL_ANG', 'HYP_DIA', 'TAP_DIA',
                     'MAX_DIA_90', 'DROP_50', 'CP_DIA25', 'CP_DIA50', 'CP_DIA75', 'CP_DIA90', 'SKL_DEPTH', 'SKL_WIDTH',
                     'RTP_COUNT']):
                print('Compute RTP skeleton')
                currT = time.time()
                rtpSkel, crownT['RTP_COUNT'], crownT['TD_MED'], crownT['TD_AVG'], crownT['MAX_DIA_90'], rtps, tips, \
                crownT['SKL_WIDTH'], crownT['SKL_DEPTH'] = rtp.getRTPSkeleton(path, skelGraph, True)
                seg.setTips(tips)
                print('RTP Skeleton computed ' + str(time.time() - currT) + 's')

            allPara[len(allPara) - 1][2] = seg.getFail()

            if ifAnyKeyIsTrue(['RDISTR_X', 'RDISTR_Y']):
                print('Compute spatial root distribution')
                currT = time.time()
                crownT['RDISTR_X'], crownT['RDISTR_Y'] = analysis.getSymmetry(rtps, rtpSkel)
                print('Symmetry computed ' + str(time.time() - currT) + 's')

            if rtpSkel != -1:
                if ifAnyKeyIsTrue(
                        ['NR_RTP_SEG_I', 'NR_RTP_SEG_II', 'ADVT_COUNT', 'BASAL_COUNT', 'ADVT_ANG', 'BASAL_ANG',
                         'HYP_DIA', 'TAP_DIA']):
                    print('searching for hypocotyl')
                    currT = time.time()
                    branchRad, nrPaths = seg.findHypocotylCluster(path, rtpSkel)
                    print('hypocotyl computed ' + str(time.time() - currT) + 's')
                    print('starting kmeans')
                    try:
                        currT = time.time()
                        c1x, c1y, c2x, c2y = analysis.plotDiaRadius(nrPaths, branchRad, path, 2)

                        print('2 clusters computed in ' + str(time.time() - currT) + 's')

                        currT = time.time()
                        segImg = seg.makeSegmentationPicture(path, rtpSkel, img, xScale, yScale, c1x, c1y, c2x, c2y)
                        imageio.imwrite(io.getFileName() + 'Seg2.png', segImg)
                        crownT['ADVT_COUNT'], crownT['BASAL_COUNT'], crownT['NR_RTP_SEG_I'], crownT['NR_RTP_SEG_II'], \
                        crownT['HYP_DIA'], crownT['TAP_DIA'] = analysis.countRootsPerSegment(c1y, c2y, c1x, c2x)
                    except:
                        c1x = None
                        c1y = None
                        c2x = None
                        c2y = None
                        pass
                    crownT['DROP_50'] = analysis.RTPsOverDepth(path, rtpSkel)
                    print('count roots per segment')
                    print('Root classes computed in ' + str(time.time() - currT) + 's')

                if ifAnyKeyIsTrue(
                        ['ADVT_ANG', 'BASAL_ANG', 'STA_RANGE', 'STA_DOM_I', 'STA_DOM_II', 'STA_25_I', 'STA_25_II',
                         'STA_50_I', 'STA_50_II', 'STA_75_I', 'STA_75_II', 'STA_90_I', 'STA_90_II', 'RTA_DOM_I',
                         'RTA_DOM_II', 'STA_MIN', 'STA_MAX', 'STA_MED', 'RTA_RANGE', 'RTA_MIN', 'RTA_MAX', 'RTA_MED']):
                    currT = time.time()
                    lat, corrBranchpts = seg.findLaterals(rtps, rtpSkel, (xScale + yScale) / 2, None)
                    print('seg.findLaterals computed in ' + str(time.time() - currT) + 's')
                    print('Compute angles at 2cm')
                    currT = time.time()
                    if c1x != None and c1y != None and c2x != None and c2y != None:
                        crownT['ADVT_ANG'], crownT['BASAL_ANG'] = analysis.anglesPerClusterAtDist(c1y, c2y, rtpSkel,
                                                                                                  path, lat,
                                                                                                  corrBranchpts,
                                                                                                  (xScale + yScale) / 2,
                                                                                                  dist=20)
                    else:
                        crownT['ADVT_ANG'] = 'nan'
                        crownT['BASAL_NG'] = 'nan'
                    print('angles at 2cm computed in ' + str(time.time() - currT) + 's')

                    if ifAnyKeyIsTrue(
                            ['STA_25_I', 'STA_25_II', 'STA_50_I', 'STA_50_II', 'STA_75_I', 'STA_75_II', 'STA_90_I',
                             'STA_90_II']):
                        try:
                            print('compute quantile angles')
                            currT = time.time()
                            a25, a50, a75, a90 = analysis.calculateAngleQuantiles(path, lat, corrBranchpts, rtpSkel)
                            print('angles computed in ' + str(time.time() - currT) + 's')
                        except:
                            a25 = ['nan']
                            a50 = ['nan']
                            a75 = ['nan']
                            a90 = ['nan']
                            print('ERROR: No quantile angles calculated')

                    if ifAnyKeyIsTrue(['RTA_RANGE', 'RTA_MIN', 'RTA_MAX', 'RTA_MED']):
                        try:
                            print('compute angles')
                            currT = time.time()
                            crownT['RTA_MED'], crownT['RTA_MIN'], crownT['RTA_MAX'], crownT[
                                'RTA_RANGE'], anglesN = analysis.calculateAngles(path, lat, corrBranchpts, rtpSkel)
                            print('RTA angle characteristics computed in ' + str(time.time() - currT) + 's')
                        except:
                            print('ERROR: No RTA angles calculated')

                    if ifAnyKeyIsTrue(['STA_RANGE', 'STA_MIN', 'STA_MAX', 'STA_MED']):
                        try:
                            print('compute STA angles')
                            currT = time.time()
                            crownT['STA_RANGE'], crownT['STA_MED'], crownT['STA_MIN'], crownT[
                                'STA_MAX'], angles = analysis.getLateralAngles(path, lat, corrBranchpts, rtpSkel)
                            print('STA angles characteristics computed in ' + str(time.time() - currT) + 's')
                        except:
                            print('ERROR: No STA angles calculated')

                    if ifAnyKeyIsTrue(['CP_DIA25', 'CP_DIA50', 'CP_DIA75', 'CP_DIA90']):
                        try:
                            print('compute diameter quantils')
                            currT = time.time()
                            crownT['CP_DIA25'], crownT['CP_DIA50'], crownT['CP_DIA75'], crownT[
                                'CP_DIA90'] = analysis.getDiameterQuantilesAlongSinglePath(path, rtpSkel)
                            print('Tap diameters computed in ' + str(time.time() - currT) + 's')
                        except:
                            print('ERROR: No quantile diameters calculated')

                    if ifAnyKeyIsTrue(['STA_DOM_I', 'STA_DOM_II']):
                        try:
                            print('compute STA dominant angles')
                            currT = time.time()
                            crownT['STA_DOM_I'], crownT['STA_DOM_II'] = analysis.findHistoPeaks(angles)
                            print('STA dominant angles computed in ' + str(time.time() - currT) + 's')
                        except:
                            print('ERROR: No dominant angles calculated (STA)')

                    if ifAnyKeyIsTrue(['STA_25_I', 'STA_25_II']):
                        try:
                            currT = time.time()
                            crownT['STA_25_I'], crownT['STA_25_II'] = analysis.findHistoPeaks(a25)
                            print('STA 25 angles computed in ' + str(time.time() - currT) + 's')
                        except:
                            print('ERROR: No dominant angles25 calculated')

                    if ifAnyKeyIsTrue(['STA_50_I', 'STA_50_II']):
                        try:
                            currT = time.time()
                            crownT['STA_50_I'], crownT['STA_50_II'] = analysis.findHistoPeaks(a50)
                            print('STA 50 angles computed in ' + str(time.time() - currT) + 's')
                        except:
                            print('ERROR: No dominant angles50 calculated')

                    if ifAnyKeyIsTrue(['STA_75_I', 'STA_75_II']):
                        try:
                            currT = time.time()
                            crownT['STA_75_I'], crownT['STA_75_II'] = analysis.findHistoPeaks(a75)
                            print('STA 75 angles computed in ' + str(time.time() - currT) + 's')
                        except:
                            print('ERROR: No dominant angles75 calculated')

                    if ifAnyKeyIsTrue(['STA_90_I', 'STA_90_II']):
                        try:
                            currT = time.time()
                            crownT['STA_90_I'], crownT['STA_90_II'] = analysis.findHistoPeaks(a90)
                            print('STA 90 angles computed in ' + str(time.time() - currT) + 's')
                        except:
                            print('ERROR: No dominant angles90 calculated')

                    if ifAnyKeyIsTrue(['RTA_DOM_I', 'RTA_DOM_II']):
                        try:
                            currT = time.time()
                            crownT['RTA_DOM_I'], crownT['RTA_DOM_II'] = analysis.findHistoPeaks(anglesN)
                            print('angles computed in ' + str(time.time() - currT) + 's')
                        except:
                            print('ERROR: No dominant RTA angles calculated')
    io.setHomePath(oldHome)
    if maxExRoot >= 1:
        rtpSkel = -1
        os.chdir(io.getHomePath())
        f = glob.glob('*.lateral.png')
        for (counter, i) in enumerate(f):
            print('processing lateral file: ' + i)

            if maxExRoot > 0:
                xScale = allPara[counter / maxExRoot][7]
                yScale = allPara[counter / maxExRoot][8]
                io.setFileName(os.path.basename(i))
            else:
                xScale = allPara[counter][7]
                yScale = allPara[counter][8]
                io.setFileName(os.path.basename(i))
                io.setidIdx(counter)

            rtp = RootTipPaths.RootTipPaths(io)

            analysis = Analysis.Analysis(io, (xScale + yScale) / 2)

            try:
                img = imageio.imread(i, as_gray=True)
            except:
                print('Image not readable')
                img = []
                pass
            if len(img) > 0:

                seg = Segmentation.Segmentation(img, io=io)
                imgL = seg.label()

                if imgL != None:
                    skel = Skeleton.Skeleton(imgL)
                    testSkel, testDia = skel.skel(imgL)
                    path, skelGraph = seg.findThickestPathLateral(testSkel, testDia, xScale, yScale)
                    if ifAnyKeyIsTrue(
                            ['LT_AVG_LEN', 'NODAL_LEN', 'LT_BRA_FRQ', 'NODAL_AVG_DIA', 'LT_AVG_ANG', 'LT_ANG_RANGE',
                             'LT_MIN_ANG', 'LT_MAX_ANG', 'LT_DIST_FIRST', 'LT_MED_DIA', 'LT_AVG_DIA']):
                        rtpSkel, _, crownT['LT_MED_DIA'], crownT['LT_AVG_DIA'], _, rtps, _, _, _ = rtp.getRTPSkeleton(
                            path, skelGraph, True)

                    if rtpSkel != -1:
                        if ifAnyKeyIsTrue(['LT_BRA_FRQ']):
                            crownT['LT_BRA_FRQ'] = analysis.getBranchingfrequencyAlongSinglePath(rtps, path)
                            crownT['NODAL_AVG_DIA'], _ = analysis.getDiametersAlongSinglePath(path, rtpSkel,
                                                                                              (xScale + yScale) / 2)
                            crownT['NODAL_LEN'] = analysis.getLengthOfPath(path)
                        if ifAnyKeyIsTrue(['LT_DIST_FIRST', 'LT_AVG_LEN', 'LT_BRA_FRQ', 'LT_ANG_RANGE', 'LT_AVG_ANG',
                                           'LT_MIN_ANG', 'LT_MAX_ANG']):
                            lat, corrBranchpts, crownT['LT_DIST_FIRST'] = seg.findLaterals(rtps, rtpSkel,
                                                                                           (xScale + yScale) / 2, path)
                            if ifAnyKeyIsTrue(['LT_AVG_LEN']):
                                crownT['LT_AVG_LEN'] = analysis.getLateralLength(lat, path, rtpSkel)
                            if ifAnyKeyIsTrue(['LT_ANG_RANGE', 'LT_AVG_ANG', 'LT_MIN_ANG', 'LT_MAX_ANG']):
                                crownT['LT_ANG_RANGE'], crownT['LT_AVG_ANG'], crownT['LT_MIN_ANG'], crownT[
                                    'LT_MAX_ANG'], _ = analysis.getLateralAngles(path, lat, corrBranchpts, rtpSkel)
            allCrown.append(crownT.copy())
    else:
        allCrown.append(crownT.copy())

    io.setHomePath(oldHome)
    # os.chdir('../')


def printHeader():
    if not os.path.exists('./options.csv') and len(sys.argv) != 13:
        print('------------------------------------------------------------')
        print('DIRT 1.1 - An automatic highthroughput root phenotyping platform')
        print('(c) 2014 Alexander Bucksch - bucksch@uga.edu')
        print('Web application by Abhiram Das - abhiram.das@gmail.com')
        print(' ')
        print('http://dirt.iplantcollaborative.org')
        print(' ')
        print('University of Georgia')
        print('------------------------------------------------------------')
        print('Program usage: python main.py (please configure the program with the otions.csv file)')
        print('<run file path> full path to file with the root image')
        print('<unique id> ID which will be a folder name in theworking directory. Integer value needed')
        print(
            '<mask threshold> multiplier for the automatically determned mask threshold. 1.0 works fine and is default. If flashlight is used, the 0.6 is a good choice.')
        print('<excised root> 1 - excised root analysis is on, 0 - excised root analysis is off')
        print('<crown root> 1 - crown root analysis is on, 0 - crown root analysis is off')
        print('<segmentation> 1 -  is on, 0 - is off')
        print('<marker diameter> a simple decimal e.g. 25.4. If 0.0 is used, then the output will have pixels as unit.')
        print('<stem reconstruction> 1 - reconstruction is turned on, 0 - reconstruction is turned off')
        print('<plots> 1 - plotting data is stored, 0 - plotting data is not stored')
        print(
            '<output format> 1 - the full trait set is put into one excel file containing empty cells for traits that were not computed, 0 - only computed files are written to the output file')
        print('<working directory> full path to folder were the result is stored')
        print('<trait file path> full path to .csv file containing the traits to be computed')
        print(' ')
        print('Example: ')
        print('/Documents/image_name.jpg 8 25.0 1 1 1 25.1 0 0 0 /Documents/image_folder/ /Documents/traits.csv')

        sys.exit()
    else:
        print('------------------------------------------------------------')
        print('DIRT 1.1 - An automatic highthroughput root phenotyping platform')
        print('(c) 2014 Alexander Bucksch - bucksch@uga.edu')
        print('Web application by Abhiram Das - abhiram.das@gmail.com')
        print(' ')
        print('http://dirt.iplantcollaborative.org')
        print(' ')
        print('University of Georgia')
        print('------------------------------------------------------------')
        print(' ')
        print('Initializing folder structure')


def main(opt=None):
    global io
    global ID
    global scale
    global allPara
    global allLat
    global allCrown
    global options
    global maxExRoot

    printHeader()

    allStart = time.time()

    if opt is None:
        options = readOptions()
    else:
        options = opt

    ID = int(options[2][1])

    try:
        scale = float(options[7][1])
    except:
        scale = 1.
    rootCrown = int(options[5][1])
    maxExRoot = int(options[4][1])
    io.__init__(options[0][1], ID=ID, plots=bool(int(options[9][1])))
    tempdir = os.path.join(options[11][1], os.path.splitext(options[1][1])[0])
    print(f"Using temp dir {tempdir}")
    init(tempdir, io)

    # Run analysis
    if int(options[6][1]) == 0:
        #io.setHomePath(os.path.join(options[11][1], str(ID)))
        print(os.getcwd())
        infile = open(os.path.join(io.getHomePath(), f"{options[1][1]}.para.sav"), 'rb')
        allPara = pickle.load(infile)
        infile.close()
        print('Saved parameters loaded')
        infile.close()

    elif int(options[6][1]) == 1:
        threadSegmentation(options[11][1], options[1][1], ID, int(options[4][1]), rootCrown, float(options[7][1]) > 0.0)
        outfile = open(os.path.join(io.getHomePath(), f"{options[1][1]}.para.sav"), 'wb')
        pickle.dump(allPara, outfile)
        outfile.close()
    else:
        print('The segmentation switch must be 0 or 1')

    if int(options[5][1]) != 0 or int(options[4][1]) != 0:
        print('Start Root Analysis')
        threadCrown(os.path.join(options[11][1], str(ID)))
        print("Exiting Root Analysis")

    compTime = int((time.time() - allStart))
    print('All done in just ' + str(compTime) + ' s!')
    print('Write output.csv file')
    r = len(allCrown)
    if r == 0: r = len(allCrown)
    for i in range(r):
        allPara[i][9] = compTime
        io.writeFile(allPara[i], allCrown[i], traitDict, int(options[10][1]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
