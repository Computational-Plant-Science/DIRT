#! /nv/hp10/adas30/bin/python
"""
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
"""
import csv
from math import ceil

from options import DIRTOptions

import imageio
import skimage
import click

from segmentation import Segmentation
from preprocessing import preprocess, PreprocessingResult
from skeletonization import skel
from analysis import Analysis
from root_tips import RootTips

import glob
import os
import time
from collections import OrderedDict

from utils import any_true, print_header, get_traits


def crown(options: DIRTOptions, preprocessing_result: PreprocessingResult):
    rtp_skel = -1
    crown_t = OrderedDict()
    crown_file_path = f"{options.input_stem}.crown.png"

    if not os.path.isfile(crown_file_path):
        raise FileNotFoundError(crown_file_path)

    print(f"Processing crown file: {crown_file_path}")

    analysis = Analysis(options, (preprocessing_result.x_scale + preprocessing_result.y_scale) / 2)
    rtp = RootTips(options)
    img = imageio.imread(crown_file_path, as_gray=True)
    segmentation = Segmentation(img)
    img_label = segmentation.label()

    crown_traits = ['AVG_DENSITY', 'WIDTH_MED', 'WIDTH_MAX', 'DIA_STM_SIMPLE', 'D10', 'D20', 'D30', 'D40', 'D50',
                    'D60', 'D70', 'D80', 'D90', 'DS10', 'DS20', 'DS30', 'DS40', 'DS50', 'DS60', 'DS70', 'DS80', 'DS90',
                    'AREA', 'ANG_TOP', 'ANG_BTM']
    if any(trait in options.traits for trait in crown_traits):
        start = time.time()
        for (trait, value) in enumerate(zip(crown_traits, list(analysis.getWidthOverHeight(
                img_label,
                preprocessing_result.x_scale,
                preprocessing_result.y_scale)))):
            crown_t[trait] = value
        print(f"Mask traits computed in {ceil(time.time() - start)} seconds")

    medial_traits = ['DIA_STM', 'TD_MED', 'TD_AVG', 'STA_RANGE', 'STA_DOM_I', 'STA_DOM_II', 'STA_25_I', 'STA_25_II',
                     'STA_50_I', 'STA_50_II', 'STA_75_I', 'STA_75_II', 'STA_90_I', 'STA_90_II', 'RTA_DOM_I',
                     'RTA_DOM_II', 'STA_MIN', 'STA_MAX', 'STA_MED', 'RTA_RANGE', 'RTA_MIN', 'RTA_MAX', 'RTA_MED',
                     'NR_RTP_SEG_I', 'NR_RTP_SEG_II', 'ADVT_COUNT', 'BASAL_COUNT', 'ADVT_ANG', 'BASAL_ANG', 'HYP_DIA',
                     'TAP_DIA', 'MAX_DIA_90', 'DROP_50', 'CP_DIA25', 'CP_DIA50', 'CP_DIA75', 'CP_DIA90', 'SKL_DEPTH',
                     'SKL_WIDTH']
    if any(trait in options.traits for trait in medial_traits):
        start = time.time()
        testSkel, testDia = skel(img_label)
        imageio.imwrite(f"{options.input_stem}.skel.png", skimage.img_as_uint(testSkel))
        print(f"Medial axis computed in {ceil(time.time() - start)} seconds")

        start = time.time()
        path, skelGraph, crown_t['DIA_STM'], skelSize = segmentation.findThickestPath(
            testSkel,
            testDia,
            preprocessing_result.x_scale,
            preprocessing_result.y_scale)
        print(f"Central path computed in {ceil(time.time() - start)} seconds")

    rtp_traits = ['TD_MED', 'TD_AVG', 'STA_RANGE', 'STA_DOM_I', 'STA_DOM_II', 'STA_25_I', 'STA_25_II', 'STA_50_I',
                  'STA_50_II', 'STA_75_I', 'STA_75_II', 'STA_90_I', 'STA_90_II', 'RTA_DOM_I', 'RTA_DOM_II',
                  'STA_MIN', 'STA_MAX', 'STA_MED', 'RTA_RANGE', 'RTA_MIN', 'RTA_MAX', 'RTA_MED', 'NR_RTP_SEG_I',
                  'NR_RTP_SEG_II', 'ADVT_COUNT', 'BASAL_COUNT', 'ADVT_ANG', 'BASAL_ANG', 'HYP_DIA', 'TAP_DIA',
                  'MAX_DIA_90', 'DROP_50', 'CP_DIA25', 'CP_DIA50', 'CP_DIA75', 'CP_DIA90', 'SKL_DEPTH', 'SKL_WIDTH',
                  'RTP_COUNT']
    if any(trait in options.traits for trait in rtp_traits):
        start = time.time()
        rtp_skel, crown_t['RTP_COUNT'], crown_t['TD_MED'], crown_t['TD_AVG'], crown_t['MAX_DIA_90'], rtps, tips, \
        crown_t['SKL_WIDTH'], crown_t['SKL_DEPTH'] = rtp.get_rtp_skeleton(path, skelGraph, True)
        segmentation.setTips(tips)
        print(f"RTP skeleton computed in {ceil(time.time() - start)} seconds")

    srd_traits = ['RDISTR_X', 'RDISTR_Y']
    if any(trait in options.traits for trait in srd_traits):
        start = time.time()
        crown_t['RDISTR_X'], crown_t['RDISTR_Y'] = analysis.get_symmetry(rtps, rtp_skel)
        print(f"Spatial root distribution computed in {ceil(time.time() - start)} seconds")

    hypocotol_traits = ['NR_RTP_SEG_I', 'NR_RTP_SEG_II', 'ADVT_COUNT', 'BASAL_COUNT', 'ADVT_ANG', 'BASAL_ANG',
                         'HYP_DIA', 'TAP_DIA']
    if rtp_skel != -1 and any(trait in options.traits for trait in hypocotol_traits):
        start = time.time()
        branchRad, nrPaths = segmentation.findHypocotylCluster(path, rtp_skel)
        print(f"Hypocotol computed in {ceil(time.time() - start)} seconds")

        try:
            start = time.time()
            c1x, c1y, c2x, c2y = analysis.plotDiaRadius(nrPaths, branchRad, path, 2)
            print(f"2 k-means clusters computed in {ceil(time.time() - start)} seconds")

            start = time.time()
            segImg = segmentation.makeSegmentationPicture(path, rtp_skel, img, preprocessing_result.x_scale, preprocessing_result.y_scale, c1x, c1y, c2x, c2y)
            imageio.imwrite(f"{options.input_stem}.seg2.png", segImg)
            crown_t['ADVT_COUNT'], crown_t['BASAL_COUNT'], crown_t['NR_RTP_SEG_I'], crown_t[
                'NR_RTP_SEG_II'], crown_t['HYP_DIA'], crown_t['TAP_DIA'] = analysis.countRootsPerSegment(c1y, c2y, c1x, c2x)
        except:
            c1x = None
            c1y = None
            c2x = None
            c2y = None
            pass
        crown_t['DROP_50'] = analysis.RTPsOverDepth(path, rtp_skel)
        print(f"Root classes computed in {ceil(time.time() - start)} seconds")

    lateral_traits = ['ADVT_ANG', 'BASAL_ANG', 'STA_RANGE', 'STA_DOM_I', 'STA_DOM_II', 'STA_25_I', 'STA_25_II',
                         'STA_50_I', 'STA_50_II', 'STA_75_I', 'STA_75_II', 'STA_90_I', 'STA_90_II', 'RTA_DOM_I',
                         'RTA_DOM_II', 'STA_MIN', 'STA_MAX', 'STA_MED', 'RTA_RANGE', 'RTA_MIN', 'RTA_MAX', 'RTA_MED']
    if rtp_skel != -1 and any(trait in options.traits for trait in lateral_traits):
        start = time.time()
        lat, corrBranchpts = segmentation.findLaterals(
            rtps,
            rtp_skel,
            (preprocessing_result.x_scale + preprocessing_result.y_scale) / 2,
            None)
        print(f"Laterals computed in {ceil(time.time() - start)} seconds")

        start = time.time()
        if c1x != None and c1y != None and c2x != None and c2y != None:
            crown_t['ADVT_ANG'], crown_t['BASAL_ANG'] = analysis.anglesPerClusterAtDist(
                c1y,
                c2y,
                rtp_skel,
                path,
                lat,
                corrBranchpts,
                (preprocessing_result.x_scale + preprocessing_result.y_scale) / 2,
                dist=20)
        else:
            crown_t['ADVT_ANG'] = 'nan'
            crown_t['BASAL_NG'] = 'nan'
        print(f"Angles at 2cm computed in {ceil(time.time() - start)} seconds")

    quantile_traits = ['STA_25_I', 'STA_25_II', 'STA_50_I', 'STA_50_II', 'STA_75_I', 'STA_75_II', 'STA_90_I', 'STA_90_II']
    if rtp_skel != -1 and any(trait in options.traits for trait in quantile_traits):
        try:
            start = time.time()
            a25, a50, a75, a90 = analysis.calculateAngleQuantiles(path, lat, corrBranchpts, rtp_skel)
            print(f"Quantile angles computed in {ceil(time.time() - start)} seconds")
        except:
            a25 = ['nan']
            a50 = ['nan']
            a75 = ['nan']
            a90 = ['nan']
            print('No quantile angles calculated')

    rta_traits = ['RTA_RANGE', 'RTA_MIN', 'RTA_MAX', 'RTA_MED']
    if rtp_skel != -1 and any(trait in options.traits for trait in rta_traits):
        start = time.time()
        crown_t['RTA_MED'], crown_t['RTA_MIN'], crown_t['RTA_MAX'], crown_t[
            'RTA_RANGE'], anglesN = analysis.calculateAngles(path, lat, corrBranchpts, rtp_skel)
        print(f"RTA angle characteristics computed in {ceil(time.time() - start)} seconds")

    sta_traits = ['STA_RANGE', 'STA_MIN', 'STA_MAX', 'STA_MED']
    if rtp_skel != -1 and any(trait in options.traits for trait in sta_traits):
        start = time.time()
        crown_t['STA_RANGE'], crown_t['STA_MED'], crown_t['STA_MIN'], crown_t[
            'STA_MAX'], angles = analysis.get_lateral_angles(path, lat, corrBranchpts, rtp_skel)
        print(f"STA angle characteristics computed in {ceil(time.time() - start)} seconds")

    diameter_traits = ['CP_DIA25', 'CP_DIA50', 'CP_DIA75', 'CP_DIA90']
    if rtp_skel != -1 and any(trait in options.traits for trait in diameter_traits):
        start = time.time()
        crown_t['CP_DIA25'], crown_t['CP_DIA50'], crown_t['CP_DIA75'], crown_t[
            'CP_DIA90'] = analysis.getDiameterQuantilesAlongSinglePath(path, rtp_skel)
        print(f"Tap diameters computed in {ceil(time.time() - start)} seconds")

    sta_dominant_traits = ['STA_DOM_I', 'STA_DOM_II']
    if rtp_skel != -1 and any(trait in options.traits for trait in sta_dominant_traits):
        start = time.time()
        crown_t['STA_DOM_I'], crown_t['STA_DOM_II'] = analysis.find_histo_peaks(angles)
        print(f"STA dominant angles computed in {ceil(time.time() - start)} seconds")

    dominant_25_traits = ['STA_25_I', 'STA_25_II']
    if rtp_skel != -1 and any(trait in options.traits for trait in dominant_25_traits):
        start = time.time()
        crown_t['STA_25_I'], crown_t['STA_25_II'] = analysis.find_histo_peaks(a25)
        print(f"STA 25 angles computed in {ceil(time.time() - start)} seconds")

    dominant_50_traits = ['STA_50_I', 'STA_50_II']
    if rtp_skel != -1 and any(trait in options.traits for trait in dominant_50_traits):
        start = time.time()
        crown_t['STA_50_I'], crown_t['STA_50_II'] = analysis.find_histo_peaks(a50)
        print(f"STA 50 angles computed in {ceil(time.time() - start)} seconds")

    dominant_75_traits = ['STA_75_I', 'STA_75_II']
    if rtp_skel != -1 and any(trait in options.traits for trait in dominant_75_traits):
        start = time.time()
        crown_t['STA_75_I'], crown_t['STA_75_II'] = analysis.find_histo_peaks(a75)
        print(f"STA 75 angles computed in {ceil(time.time() - start)} seconds")

    dominant_90_traits = ['STA_90_I', 'STA_90_II']
    if rtp_skel != -1 and any(trait in options.traits for trait in dominant_90_traits):
        start = time.time()
        crown_t['STA_90_I'], crown_t['STA_90_II'] = analysis.find_histo_peaks(a90)
        print(f"STA 90 angles computed in {ceil(time.time() - start)} seconds")

    rta_dominant_traits = ['RTA_DOM_I', 'RTA_DOM_II']
    if rtp_skel != -1 and any(trait in options.traits for trait in rta_dominant_traits):
        start = time.time()
        crown_t['RTA_DOM_I'], crown_t['RTA_DOM_II'] = analysis.find_histo_peaks(anglesN)
        print(f"RTA dominant angles computed in {ceil(time.time() - start)} seconds")

    if options.excised_roots == 0:
        return crown_t

    rtp_skel = -1
    lateral_file_path = f"{options.input_stem}.crown.png"
    rtp = RootTips(options)
    analyze = Analysis(options, (preprocessing_result.x_scale + preprocessing_result.y_scale) / 2)
    img = imageio.imread(lateral_file_path, as_gray=True)
    segmentation = Segmentation(img)
    image_label = segmentation.label()

    if image_label is None:
        raise ValueError(f"No image label found, skipping lateral file processing")
    print(f"Processing lateral file: {lateral_file_path}")

    testSkel, testDia = skel(image_label)
    path, skelGraph = segmentation.findThickestPathLateral(
        testSkel,
        testDia,
        preprocessing_result.x_scale,
        preprocessing_result.y_scale)

    rtp_skeleton_traits = ['LT_AVG_LEN', 'NODAL_LEN', 'LT_BRA_FRQ', 'NODAL_AVG_DIA', 'LT_AVG_ANG', 'LT_ANG_RANGE',
             'LT_MIN_ANG', 'LT_MAX_ANG', 'LT_DIST_FIRST', 'LT_MED_DIA', 'LT_AVG_DIA']
    if any(trait in options.traits for trait in rtp_skeleton_traits):
        start = time.time()
        rtp_skel, _, crown_t['LT_MED_DIA'], crown_t[
            'LT_AVG_DIA'], _, rtps, _, _, _ = rtp.get_rtp_skeleton(
            path, skelGraph, True)
        print(f"RTP skeleton computed in {ceil(time.time() - start)} seconds")

    branching_freq_traits = ['LT_BRA_FRQ']
    if rtp_skel != -1 and any(trait in options.traits for trait in branching_freq_traits):
        start = time.time()
        crown_t['LT_BRA_FRQ'] = analyze.get_branching_frequency_along_single_path(rtps, path)
        crown_t['NODAL_AVG_DIA'], _ = analyze.getDiametersAlongSinglePath(
            path,
            rtp_skel,
            (preprocessing_result.x_scale + preprocessing_result.y_scale) / 2)
        crown_t['NODAL_LEN'] = analyze.getLengthOfPath(path)
        print(f"Branching frequency computed in {ceil(time.time() - start)} seconds")

    rtp_lateral_traits = ['LT_DIST_FIRST', 'LT_AVG_LEN', 'LT_BRA_FRQ', 'LT_ANG_RANGE', 'LT_AVG_ANG',
                     'LT_MIN_ANG', 'LT_MAX_ANG']
    if rtp_skel != 1 and any(trait in options.traits for trait in rtp_lateral_traits):
        start = time.time()
        lat, corrBranchpts, crown_t['LT_DIST_FIRST'] = segmentation.findLaterals(
            rtps,
            rtp_skel,
            (preprocessing_result.x_scale + preprocessing_result.y_scale) / 2,
            path)
        print(f"Laterals computed in {ceil(time.time() - start)} seconds")

    lateral_length_traits = ['LT_AVG_LEN']
    if rtp_skel != 1 and any(trait in options.traits for trait in lateral_length_traits):
        start = time.time()
        crown_t['LT_AVG_LEN'] = analyze.get_lateral_length(lat, path, rtp_skel)
        print(f"Lateral length computed in {ceil(time.time() - start)} seconds")

    lateral_angles_traits = ['LT_ANG_RANGE', 'LT_AVG_ANG', 'LT_MIN_ANG', 'LT_MAX_ANG']
    if rtp_skel != 1 and any(trait in options.traits for trait in lateral_angles_traits):
        start = time.time()
        crown_t['LT_ANG_RANGE'], crown_t['LT_AVG_ANG'], crown_t['LT_MIN_ANG'], crown_t[
            'LT_MAX_ANG'], _ = analyze.get_lateral_angles(path, lat, corrBranchpts, rtp_skel)
        print(f"Lateral angles computed in {ceil(time.time() - start)} seconds")


columns = [
    'Image',
    'Failed',
    'Tag',
    'Circle Ratio',
    'X Pixel',
    'Y Pixel',
    'X Scale',
    'Y Scale',
    'Duration',
    'Skeleton Vertices'
]


@click.command()
@click.argument('input_path')
@click.option('--trait_file_path', required=True, type=str)
@click.option('--excised_roots', required=False, type=int, default=0)
@click.option('--mask_threshold', required=False, type=float, default=10.0)
@click.option('--marker_diameter', required=False, type=float, default=25.4)
@click.option('--crown_root', required=False, type=bool, default=True)
@click.option('--segmentation', required=False, type=bool, default=True)
@click.option('--stem_reconstruction', required=False, type=bool, default=False)
@click.option('--plot', required=False, type=bool, default=True)
def run(input_path,
        trait_file_path,
        excised_roots,
        mask_threshold,
        marker_diameter,
        crown_root,
        segmentation,
        stem_reconstruction,
        plot):
    print_header()

    start = time.time()
    traits = get_traits(trait_file_path)
    options = DIRTOptions(
        input_path=input_path,
        trait_file_path=trait_file_path,
        excised_roots=excised_roots,
        mask_threshold=mask_threshold,
        marker_diameter=marker_diameter,
        crown_root=crown_root,
        segmentation=segmentation,
        stem_reconstruction=stem_reconstruction,
        plot=plot,
        traits=traits)

    if options.segmentation:
        preprocessing_result = preprocess(options)
        print(preprocessing_result)
        if options.excised_roots or options.crown_root:
            crown_result = crown(options, preprocessing_result)
            print(crown_result)

    with open('output.csv', 'a' if os.path.isfile('output.csv') else 'w', newline='') as output_file:
        writer = csv.writer(output_file, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(columns + traits)
        # TODO write output rows

    print(f"All done in just {ceil((time.time() - start))} seconds!")


if __name__ == '__main__':
    run()