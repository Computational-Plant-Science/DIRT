#! /nv/hp10/adas30/bin/python
'''
----------------------------------------------------------------------------------------------------
DIRT - An automatic highthroughput root phenotyping platform
Web interface by Abhiram Das - adas30@biology.gatech.edu

http://www.dirt.biology.gatech.edu

Georgia Institute of Technology

The software is written in:
- python 2.7 (https://www.python.org)

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
- adaptive threshholding in Masking.py (http://scikit-image.org)

The software uses modified code from Kyle Fox:
    - fixOrientation.py: https://github.com/kylefox/python-image-orientation-patch


Please cite the DIRT Paper if you use the code for your scientific project.

Bucksch et al., 2014 "Image-based high-throughput field phenotyping of crop roots", Plant Physiology


----------------------------------------------------------------------------------------------------
Author: Alexander Bucksch
School of Biology and Interactive computing
Georgia Institute of Technology

Mail: bucksch@gatech.edu
Web: http://www.bucksch.nl
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

'''
# internal library imports
'''
import IO
import Segmentation
import Preprocessing
import Skeleton
import Analysis
import RootTipPaths
import time
from fixImageOrientation import *

'''
# standard python imports
'''
import os
import pickle
import csv
import sys

'''
#global defs
'''
allLat=[]
allPara=[]
allCrown=[]
f=[]
imgID=None
def init(fpath,io):
    oldpath = os.getcwd()
    io.setHomePath(fpath)
    if os.path.exists(fpath)==False: 
        os.mkdir(fpath)
    os.chdir(fpath)
    print os.getcwd()
    io.setServerPath(os.getcwd())
    if os.path.exists('./tmp/') == False:
        os.mkdir('./tmp/')
    if os.path.exists('./Mask/') == False:
        os.mkdir('./Mask/')
    if os.path.exists('./Lateral/') == False:
        os.mkdir('./Lateral/')
        if os.path.exists('./Lateral/Plots/') == False:
            os.mkdir('./Lateral/Plots/')
        if os.path.exists('./Lateral/Result/') == False:
            os.mkdir('./Lateral/Result/')
    else:
        if os.path.exists('./Lateral/Plots/') == False:
            os.mkdir('./Lateral/Plots/')
        if os.path.exists('./Lateral/Result/') == False:
            os.mkdir('./Lateral/Result/')
    
    if os.path.exists('./Crown/') == False:
        os.mkdir('./Crown/')
        if os.path.exists('./Crown/Plots/') == False:
            os.mkdir('./Crown/Plots/')
        if os.path.exists('./Crown/Result/') == False:
            os.mkdir('./Crown/Result/')
    else:
        if os.path.exists('./Crown/Plots/') == False:
            os.mkdir('./Crown/Plots/')
        if os.path.exists('./Crown/Result/') == False:
            os.mkdir('./Crown/Result/')
    
    os.chdir(oldpath)
    
def readOptions():
    options=[]
    if len(sys.argv)==10:
        options.append([0,os.path.dirname(sys.argv[1])+'/'])
        options.append([0,os.path.basename(sys.argv[1])])
        options.append([0,sys.argv[2]])
        options.append([0,sys.argv[3]])
        options.append([0,sys.argv[4]])
        options.append([0,sys.argv[5]])
        options.append([0,sys.argv[6]])
        options.append([0,sys.argv[7]])
        options.append([0,sys.argv[8]])
        options.append([0,sys.argv[9]])
    else:
        with open('./options.csv','U') as csvfile: 
            filedata= csv.reader(csvfile)
            for i in filedata:
                options.append(i)
        
        
    return options

def threadSegmentation(filepath,imgFile,imgID,maxExRoot,marker):
    
    io.setFileName(imgFile)
    io.setidIdx(imgID)
    prep=Preprocessing.Preprocessing(io)
    print 'segmenting file: '+imgFile +'\n'

    if os.path.isfile(options[0][1]+imgFile):
        # fix orientation of the image in tiff and Jpg files
        fix_orientation(options[0][1]+imgFile, save_over=True)
        img= scipy.misc.imread(options[0][1]+imgFile,flatten=True)
            
    else:
        print 'Image not readable'
        img=[]

    if len(img)>0: 
        currT=time.time()       
        Failed,tagExtract,circleRatio, circleWidth, circleHeight = prep.prepocess(img,scale=float(options[3][1]),nrExRoot=maxExRoot,marker=marker)
        print 'Segmentation finished in '+str(time.time()-currT)+'s'
        if Failed == False:
            xScale=scale/float(circleWidth)
            yScale=scale/float(circleHeight)
            if xScale<=0.0: xScale=1.
            if yScale<=0.0: yScale=1.
            para=[int(imgID),io.getFileName(),Failed,tagExtract,circleRatio, circleWidth, circleHeight,xScale,yScale,-1,-1]
            if maxExRoot>1:
                for _ in range(maxExRoot):
                    allPara.append(para)
            else: allPara.append(para)
        else:
            xScale=scale/float(1.0)
            yScale=scale/float(1.0)
            circleRatio=1.0
            circleWidth=1.0
            circleHeight=1.0
            para=[int(imgID),io.getFileName(),Failed,tagExtract,circleRatio, circleWidth, circleHeight,xScale,yScale,-1,-1]
            if maxExRoot>1:
                for _ in range(maxExRoot):
                    allPara.append(para)
            else: allPara.append(para)

def threadCrown(filepath):
    imgL=[]
    tipdiameter=float(options[8][1])
    print io.getHomePath()
    os.chdir(io.getHomePath())
    io.setHomePath('./Crown/')
    f=io.scanDir()
    for (counter,i) in enumerate(f):
        io.setFileName(os.path.basename(i))
        io.setidIdx(imgID)
        
        print 'processing Crown file: '+i
        xScale=allPara[counter][7]
        yScale=allPara[counter][8]
        analysis=Analysis.Analysis(io,(xScale+yScale)/2)
        rtp=RootTipPaths.RootTipPaths(io,tp=tipdiameter)
        rtp.setTipDiaFilter(tipdiameter*(xScale+yScale)/2)
        crownT=[]
        
        try:
            img=scipy.misc.imread(i,flatten=True)
        except:
            print 'Image not readable'
            img=-1
        if len(img)>0:
            seg=Segmentation.Segmentation(img,io)
            imgL=seg.label()
            print 'compute root profile'
            currT=time.time()
            rootDensity,medianWidth,maxWidth,D,DS,_,_,_,_=analysis.getWidthOverHeight(imgL,xScale,yScale)
            print 'Mask traits computed '+str(time.time()-currT)+'s'
            currT=time.time()
            skel=Skeleton.Skeleton(imgL)
            testSkel,testDia=skel.skel(imgL)
            print 'Medial axis computed '+str(time.time()-currT)+'s'
            currT=time.time()
            path,skelGraph,stemDia,skelSize=seg.findThickestPath(testSkel,testDia,xScale,yScale)
            allPara[counter][10]=skelSize
            print 'Central path computed '+str(time.time()-currT)+'s'
            print 'compute rtp skeleton'
            currT=time.time()
            rtpSkel,nrOfRTP, medianTipDiameter,meanDiameter,dia90, _, rtps, tips, _, _ =rtp.getRTPSkeleton(path,skelGraph,True)
            allPara[len(allPara)-1][2]=seg.getFail()
            seg.setTips(tips)
            print 'RTP Skeleton computed '+str(time.time()-currT)+'s'
            print 'compute symmetry'
            currT=time.time()
            vecSym=analysis.getSymmetry(rtps,rtpSkel)
            print 'Symmetry computed '+str(time.time()-currT)+'s'
            
            if rtpSkel!=-1:

                lat,corrBranchpts,_=seg.findLaterals(rtps, rtpSkel,(xScale+yScale)/2)
                    
                try:
                    print 'compute quantile angles'
                    currT=time.time()
                    a25,a50,a75,a90=analysis.calculateAngleQuantiles(path,lat,corrBranchpts,rtpSkel)
                    print 'angles computed in '+str(time.time()-currT)+'s'
                except:
                    raise
                    a25=['nan']
                    a50=['nan']
                    a75=['nan']
                    a90=['nan']

                    print 'ERROR: No quantile angles calculated'
                
                try:
                    print 'compute angles'
                    currT=time.time()
                    angRangeN,avgAngleN,minangleN,maxAngleN,anglesN=analysis.calculateAngles(path,lat,corrBranchpts,rtpSkel)
                    print 'angles computed in '+str(time.time()-currT)+'s'
                except:
                    avgAngleN='nan'
                    minangleN='nan'
                    maxAngleN='nan'
                    angRangeN='nan'
                    anglesN='nan'
                    print 'ERROR: No angles calculated'
                    
                    
                try:
                    print 'compute RTA angles'
                    currT=time.time()
                    angRange,avgAngle,minangle,maxAngle,angles=analysis.getLateralAngles(path,lat,corrBranchpts,rtpSkel)
                    print 'angles computed in '+str(time.time()-currT)+'s'
                except:
                    raise
                    avgAngle='nan'
                    minangle='nan'
                    maxAngle='nan'
                    angRange='nan'
                    angles='nan'
                    print 'ERROR: No RTA angles calculated'
                try:

                    print 'compute diameter quantils'
                    currT=time.time()
                    d25,d50,d75,d90=analysis.getDiameterQuantilesAlongSinglePath(path,rtpSkel)
                    print 'diameters computed in '+str(time.time()-currT)+'s'
                except:
                    d25='nan'
                    d50='nan'
                    d75='nan'
                    d90='nan'
                    print 'ERROR: No quantile angles calculated'
                    raise

                
                try:
                    print 'compute dominant angles'
                    currT=time.time()
                    ang1,ang2=analysis.findHistoPeaks(angles)
                    print 'angles computed in '+str(time.time()-currT)+'s'
                except:
                    ang1='nan'
                    ang2='nan'
                    print 'ERROR: No dominant angles calculated'
                try:
                    currT=time.time() 
                    ang25_1,ang25_2=analysis.findHistoPeaks(a25)
                    print 'angles computed in '+str(time.time()-currT)+'s'
                except:
                    ang25_1='nan'
                    ang25_2='nan'
                    print 'ERROR: No dominant angles25 calculated'
                try:
                    currT=time.time()
                    ang50_1,ang50_2=analysis.findHistoPeaks(a50)
                    print 'angles computed in '+str(time.time()-currT)+'s'
                except:
                    ang50_1='nan'
                    ang50_2='nan'
                    print 'ERROR: No dominant angles50 calculated'
                try:
                    currT=time.time()     
                    ang75_1,ang75_2=analysis.findHistoPeaks(a75)
                    print 'angles computed in '+str(time.time()-currT)+'s'
                except:
                    ang75_1='nan'
                    ang75_2='nan'
                    print 'ERROR: No dominant angles75 calculated'
                try:
                    currT=time.time()    
                    ang90_1,ang90_2=analysis.findHistoPeaks(a90)
                    print 'angles computed in '+str(time.time()-currT)+'s'
                except: 
                    ang90_1='nan'
                    ang90_2='nan'
                    print 'ERROR: No dominant angles90 calculated'
                
                try:   
                    currT=time.time() 
                    angN_1,angN_2=analysis.findHistoPeaks(anglesN)
                    print 'angles computed in '+str(time.time()-currT)+'s'
                except: 
                    angN_1='nan'
                    angN_2='nan'
                    print 'ERROR: No dominant angles90 calculated'
                
                crownT=[stemDia,rootDensity,angRange,ang1,ang2,ang25_1,ang25_2,ang50_1,ang50_2,ang75_1,ang75_2,ang90_1,ang90_2,angN_1,angN_2,minangle,maxAngle,avgAngle,angRangeN,avgAngleN,minangleN,maxAngleN,nrOfRTP,medianTipDiameter,meanDiameter,dia90,medianWidth,maxWidth,D[0],D[1],D[2],D[3],D[4],D[5],D[6],D[7],D[8],DS[0],DS[1],DS[2],DS[3],DS[4],DS[5],DS[6],DS[7],DS[8],vecSym[0],vecSym[1],d25,d50,d75,d90]

            else:
                crownT=[stemDia,rootDensity,'nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan',medianWidth,maxWidth,D[0],D[1],D[2],D[3],D[4],D[5],D[6],D[7],D[8],DS[0],DS[1],DS[2],DS[3],DS[4],DS[5],DS[6],DS[7],DS[8],vecSym[0],vecSym[1],d25,d50,d75,d90]

            if maxExRoot > 1:
                for i in range(maxExRoot):
                    allCrown.append(crownT)
            else:
                allCrown.append(crownT)    
            if options[4][1] == '0':
                lateralT=['nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan']
                allLat.append(lateralT)
    io.setHomePath(os.getcwd())

def threadLateral(filepath):
    tipdiameter=0.
    os.chdir(io.getHomePath())
    io.setHomePath('./Lateral/')
    f=io.scanDir()
    for (counter,i) in enumerate(f):
        print 'processing lateral file: '+i
        if maxExRoot>0:
            xScale=allPara[counter/maxExRoot][7]
            yScale=allPara[counter/maxExRoot][8]
            io.setFileName(os.path.basename(i))
        else:
            xScale=allPara[counter][7]
            yScale=allPara[counter][8] 
            io.setFileName(os.path.basename(i))
            io.setidIdx(counter)
    
        rtp=RootTipPaths.RootTipPaths(io,tipdiameter)
        rtp.setTipDiaFilter(tipdiameter)
        
        analysis=Analysis.Analysis(io,(xScale+yScale)/2)
        
        lateralT=[]

        try:
            img=scipy.misc.imread(i,flatten=True)
        except:
            print 'Image not readable'
            img=[]
            pass
        if len(img)>0:

            seg=Segmentation.Segmentation(img,io=io)
            imgL=seg.label()

            if imgL!=None:
                skel=Skeleton.Skeleton(imgL)
                testSkel,testDia=skel.skel(imgL)
                path,skelGraph=seg.findThickestPathLateral(testSkel,testDia,xScale,yScale)
                rtpSkel,_,medianD,meanD,_,_,rtps,_,_,_=rtp.getRTPSkeleton(path,skelGraph,True)
                
                if rtpSkel!=-1:
                    lBranchFreq=analysis.getBranchingfrequencyAlongSinglePath(rtps,path)
                    avgLatDiameter,slope=analysis.getDiametersAlongSinglePath(path,rtpSkel,(xScale+yScale)/2)
                    lengthNodalRoot=analysis.getLengthOfPath(path)
                    lat,corrBranchpts,distToFirst=seg.findLaterals(rtps, rtpSkel,(xScale+yScale)/2)
                    avgLLength=analysis.getLateralLength(lat,path,rtpSkel)
                    angRange,avgAngle,minangle,maxAngle,_=analysis.getLateralAngles(path,lat,corrBranchpts,rtpSkel)
                    lateralT=[avgLLength*((xScale+yScale)/2),float(lengthNodalRoot)*((xScale+yScale)/2),lBranchFreq,avgLatDiameter,slope,avgAngle,angRange,minangle,maxAngle,float(distToFirst)*((xScale+yScale)/2),medianD,meanD]
                else:
                    lateralT=['nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan']
            else:
                lateralT=['nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan']
            allLat.append(lateralT)
            if options[5][1] == '0':
                crownT=['nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan','nan']
                allCrown.append(crownT)

    io.setHomePath(os.getcwd())   

def printHeader():
    if os.path.exists('./options.csv')==False and len(sys.argv)!=10:
        print '------------------------------------------------------------' 
        print 'DIRT - An automatic highthroughput root phenotyping platform'
        print '(c) 2014 Alexander Bucksch - bucksch@gatech.edu'
        print 'Web application by Abhiram Das - adas30@biology.gatech.edu'
        print ' '
        print 'http://www.dirt.biology.gatech.edu'
        print ' '
        print 'Georgia Institute of Technology'
        print '------------------------------------------------------------'
        print 'Program usage: python main.py (please configure the program with the otions.csv file)'
        print '<run file path> full path to file with the root image'
        print '<unique id> ID which will be a folder name in theworking directory. Integer value needed'
        print '<mask threshold> multiplier for the automatically determned mask threshold. 1.0 works fine and is default. If flashlight is used, the 0.6 is a good choice.' 
        print '<excised roots> number of roots placed at the right of the root crown, 0 - excised root analysis is off' 
        print '<crown root> 1 - crown root analysis is on, 0 - crown root analysis is off' 
        print '<segmentation> 1 -  is on, 0 - is off' 
        print '<marker diameter> a simple decimal e.g. 25.4. If 0.0 is used, then the output will have pixels as unit.'
        print '<tip diameter filter> not active anymore, but can be used to consider only paths to tips of a certain size. We suggest to use 0'
        print '<working directory> full path to folder were the result is stored'
        print ' '
        print 'Example: '
        print '/Users/image_folder/image_name.jpg 8 25.0 1 1 1 25.1 0 /Users/output_folder/'
        
        sys.exit()
    else:
        print '------------------------------------------------------------' 
        print 'DIRT - An automatic highthroughput root phenotyping platform'
        print '(c) 2014 Alexander Bucksch - bucksch@gatech.edu'
        print 'Web application by Abhiram Das - adas30@biology.gatech.edu'
        print ' '
        print 'http://www.dirt.biology.gatech.edu'
        print ' '
        print 'Georgia Institute of Technology'
        print '------------------------------------------------------------'
        print ' '
        print 'Initializing folder structure'  
       
if __name__ == '__main__':
    
    allStart=time.time()
    
    printHeader()

    #defs
    options=readOptions()
    f=options[0][1]+options[1][1]
    ID=int(options[2][1])
    
    try: maskScaleThresh=float(options[3][1])
    except: maskScaleThresh=10.
    try: scale = float(options[7][1])
    except: scale =1.
    xScale=1
    yScale=1
    maxExRoot=int(options[4][1])
    marker=True
    if float(options[7][1])==0: marker=False 
    io=IO.IO(options[0][1],ID=ID)
    init(options[9][1]+str(ID)+'/',io)
    
    #Run analysis
    if int(options[6][1]) == 0:
        tipdiameter=0.
        io.setHomePath(options[9][1]+str(ID)+'/')
        print os.getcwd()
        infile=open(io.getHomePath()+'/tmp/para.sav','rb')
        allPara=pickle.load(infile)
        infile.close()
        print 'Saved paremeters loaded'
        infile.close()
        
    elif int(options[6][1]) == 1:
        threadSegmentation(options[9][1],options[1][1],ID,int(options[4][1]),float(options[7][1])>0.0)
        outfile=open(io.getHomePath()+'/tmp/para.sav','wb')
        pickle.dump(allPara,outfile)
        outfile.close()
    else: print'The segmentation switch must be 0 or 1'
    
    if int(options[5][1]) == 1 and int(options[4][1]) == 0: 
        
        print 'Start Crown Analysis'
        threadCrown(options[9][1]+str(ID)+'/')
        print "Exiting Crown Analysis"
        
    if int(options[4][1]) >= 1 and int(options[5][1]) == 0:
        print 'Start Lateral Analysis'
        threadLateral(options[9][1]+str(ID)+'/')
        print "Exiting Lateral Analysis"
        
    if int(options[5][1]) == 1 and int(options[4][1]) >= 1: 
        print 'Start Crown Analysis'
        threadCrown(options[9][1]+str(ID)+'/')
        print 'Start Lateral Analysis'
        threadLateral(options[9][1]+str(ID)+'/')
        print "Exiting LAteral Analysis"
    
    compTime=int((time.time()-allStart))
    print 'All done in just '+str(compTime)+' s!'  
    print 'Write output.csv file'
    r=len(allLat)
    if r==0: r=len(allCrown)
    for i in range(r):
        allPara[i][9]=compTime
        io.writeFile(allPara[i], allCrown[i], allLat[i])
  
