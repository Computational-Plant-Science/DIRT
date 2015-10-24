'''
Preprocessing.py

This module contains all imaging operations to create input for analysis, 
such as detecting the root in the image create and filter the binary image etc. 

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

'''
# external library imports
'''
import numpy as np
import scipy.misc
import scipy.ndimage

'''
# internal library imports
'''
import Segmentation
import Masking
import DirtOcr as ocr

'''
# standard python imports
'''
import os

class Preprocessing(object):
    '''
    classdocs
    '''


    def __init__(self,io):
        '''
        Constructor
        '''
        self.__io=io
        self.__labelHist=[]
        self.__id=io.getID()
        self.__currentIdx=io.getCurrentID()
        self.__compsX=[]
        self.__compsY=[]
        self.__h=0
        self.__w=0
        self.__tagCrop=10
        
    def prepocess(self,img,rootCrown,scale=1.0,nrExRoot=1, marker=True, stemCorrection=False):
        print 'starting to segment'
        rIdx=-1
        self.__io.setServerPath('./')
        circleIdx= circleRatio= circleWidth= circleHeight= imgCircle = 0
        Failed=False
        orig=img.copy()
        mask=Masking.Masking(scale=scale)
        imgGrey = img.astype(np.uint8)
        print 'make mask'
        imgBinary=mask.calculateMask(imgGrey)
        print 'saving binary mask'
        scipy.misc.imsave(self.__io.getHomePath()+'/Mask/' + self.__io.getFileName()+'.png', imgBinary)
        pathold=os.getcwd()
        os.chdir(self.__io.getHomePath())
        self.__io.writeServerFile('dirt_out.csv',self.__io.getHomePath()+'/Mask/'+self.__io.getFileName()+'.png,' +str(self.__io.getID())+',0')
        os.chdir(pathold)
        imgLabel=self.calculateLabelHist(imgBinary)

        if marker== True: 
            print 'Marker is True'
            circleIdx, circleRatio, circleWidth, circleHeight, imgCircle =self.findCircle(imgLabel.copy())
        else: 
            print 'Marker is False'
            circleIdx, circleRatio, circleWidth, circleHeight, imgCircle = -1, 1, 1, 1, None
        
        rectIdx, _, _, _,imgTag, tagText =self.findTag(imgLabel.copy() , imgBinary, orig, rect_ratio=5.)
       
        if rectIdx >=0:
            print 'tagIdx'+str(rectIdx)
            try: self.__labelHist[rectIdx] = 0
            except: pass

        if circleIdx >=0 and marker==True:
            try: self.__labelHist[circleIdx] = 0
            except: pass

        '''
        These two functions belong together and have to be called right after each other. I know, that is bad.
        '''
        if rootCrown==True:
            rIdx,rIdxList,crownMin,crownMax,crownBottom,crownTop=self.findRoot(imgLabel.copy()) 
            if stemCorrection== True: 
                print 'Stem reconstruction is active '
                imgRoot=self.correctForStem(imgLabel.copy(), [circleIdx,rectIdx,rIdx], crownMin, crownMax, crownBottom, crownTop, rIdx, rIdxList)
            else:
                print 'No stem reconstruction active' 
                imgReturn=np.zeros_like(imgLabel)
                imgReturn[rIdxList]=1
                imgRoot=imgReturn[crownMax:crownMin,crownBottom:crownTop]
        
        if nrExRoot >1 and rootCrown==True:

            for i in range(nrExRoot): 
                exRIdx,imgExRoot,centerPtx,centerPty=self.findExcisedRoot(imgLabel.copy(),[circleIdx,rectIdx,rIdx],crownMin,crownMax)
                if exRIdx != -1:
                    print 'found excised root '+str(i)
                    try: 
                        scipy.misc.imsave(self.__io.getHomePath()+'/Lateral/' + self.__io.getFileName()+'_'+str(centerPtx)+'_'+str(centerPty)+'.png', imgExRoot)
                        print 'excised root '+str(i)+'saved'
                    except:
                        print 'NOT SAVED !!!'
                        raise
                    try: 
                        pathold=os.getcwd()
                        os.chdir(self.__io.getHomePath())
                        self.__io.writeServerFile('dirt_out.csv',self.__io.getHomePath()+'/Lateral/'+self.__io.getFileName()+'_'+str(centerPtx)+'_'+str(centerPty)+'.png,' +str(self.__io.getID())+',0')
                        print 'excised root '+str(i)+'saved Server'
                        os.chdir(pathold)
                    except: 
                        print 'NOT SAVED !!!!'
                        raise
        elif nrExRoot ==1 and rootCrown==True: 
            exRIdx,imgExRoot,centerPtx,centerPty=self.findExcisedRoot(imgLabel.copy(),[circleIdx,rectIdx,rIdx],crownMin,crownMax)
            if exRIdx != -1:
                print 'found the excised root '
                try: 
                    pathold=os.getcwd()
                    os.chdir(self.__io.getHomePath())
                    scipy.misc.imsave(self.__io.getHomePath()+'/Lateral/' + self.__io.getFileName()+'_'+str(centerPtx)+'_'+str(centerPty)+'.png', imgExRoot)
                    print 'excised root saved' 
                    self.__io.writeServerFile('dirt_out.csv',self.__io.getHomePath()+'/Lateral/'+self.__io.getFileName()+'_'+str(centerPtx)+'_'+str(centerPty)+'.png,' +str(self.__io.getID())+',0')
                    print 'excised root saved Server'
                    os.chdir(pathold)
                except: print 'NOT SAVED !!!!'
        elif nrExRoot ==1 and rootCrown==False:
            exRIdx,imgExRoot,centerPtx,centerPty=self.findExcisedRoot(imgLabel.copy(),[circleIdx,rectIdx],0,1)
            if exRIdx != -1:
                print 'found the excised root '
                rIdx=-1
                try: 
                    pathold=os.getcwd()
                    os.chdir(self.__io.getHomePath())
                    scipy.misc.imsave(self.__io.getHomePath()+'/Lateral/' + self.__io.getFileName()+'_'+str(centerPtx)+'_'+str(centerPty)+'.png', imgExRoot)
                    print 'excised root saved' 
                    self.__io.writeServerFile('dirt_out.csv',self.__io.getHomePath()+'/Lateral/'+self.__io.getFileName()+'_'+str(centerPtx)+'_'+str(centerPty)+'.png,' +str(self.__io.getID())+',0')
                    print 'excised root saved Server'
                    os.chdir(pathold)
                except: print 'NOT SAVED !!!!'
            
        
        if marker==True:
            scipy.misc.imsave(self.__io.getHomePath()+'/Mask/' + self.__io.getFileName()+'Circle.png', imgCircle)
            scipy.misc.imsave(self.__io.getHomePath()+'/Mask/' + self.__io.getFileName()+'Tag.png', imgTag)
        #pathold=os.getcwd()
        #os.chdir(self.__io.getHomePath())
        
        if marker==True: 
            self.__io.writeServerFile('dirt_out.csv',self.__io.getHomePath()+'/Mask/'+self.__io.getFileName()+'Circle.png,' +str(self.__io.getID())+',0')
            
        #os.chdir(pathold)
        if rIdx != -1:
            '''
            If image is usable, then it gets segmented and copied. Otherwise we ignore it
            '''
            try:
                print 'root image to be saved'
                scipy.misc.imsave(self.__io.getHomePath()+'/Crown/' + self.__io.getFileName()+'.png', imgRoot)
            except: 
                print 'CROWN NOT SAVED'
                raise
            try:
                pathold=os.getcwd()
                os.chdir(self.__io.getHomePath()) 
                self.__io.writeServerFile('dirt_out.csv',self.__io.getHomePath()+'/Crown/'+self.__io.getFileName()+'.png,' +str(self.__io.getID())+',0')
                os.chdir(pathold)
            except: print 'MASK NOT WRITTEN TO SERVER FILE'
        else: Failed=True
        print "old path: "+pathold
        return  Failed,tagText,circleRatio, circleWidth, circleHeight
        
    def calculateLabelHist(self,imgBinary):
        seg=Segmentation.Segmentation(imgBinary,io=self.__io)    
        labeled,_=seg.labelAll()
        x, y = np.shape(labeled)
        val = labeled.flatten()

        histo, _ = np.histogram(val, bins=np.max(labeled) + 1)
        '''
        TEST: background can have less pixels than foreground if no markers are in the image
        '''
        if len(np.unique(labeled))==2:
            nrOfwhitePx=len(np.where(imgBinary==255)[1])
            comp1 = np.max(histo)
            if comp1==nrOfwhitePx:
                comp1 = np.min(histo)
            idx1 = list(histo).index(comp1)
            histo[idx1]=0
        else:
            comp1 = np.max(histo)
            idx1 = list(histo).index(comp1)
            histo[idx1]=0
            
        for i in range(len(histo)):
            if histo[i]<100:
                histo[i]=0
        self.__labelHist = histo
        
        self.__compsX = [[] for i in range(len(self.__labelHist))]
        self.__compsY = [[] for i in range(len(self.__labelHist))]

        self.__h, self.__w = np.shape(labeled)
        for i in range(self.__w):
            for j in range(self.__h):
                self.__compsX[labeled[j][i]].append(i)
                self.__compsY[labeled[j][i]].append(j)
        return labeled
                
    def findCircle(self, labeled):
        print 'searching circle'
        ratio = []
        w,h=np.shape(labeled)
        for i in range(len(self.__compsX)):
            if self.__labelHist[i] > 0:
                xMin = np.min(self.__compsX[i])
                xMax = np.max(self.__compsX[i])
                yMin = np.min(self.__compsY[i])
                yMax = np.max(self.__compsY[i])
                nonZ=len(self.__compsX[i])
                allPx=(xMax-xMin)*(yMax-yMin)
                ratioSqualetoCircleRatio=float(nonZ)/float(allPx)
                '''
                compensates for small noisy components and small excised roots
                '''
                if float(nonZ)/float(w*h)>0.001:
                    '''
                    the inscribed circle of a bounding box fills exactly 78.64 percent -> we allow 8.64 percent variation due to noise
                    '''
                    tagRatio=(float(xMax) - float(xMin)) / (float(yMax) - float(yMin))
                    print 'Circle Ratio before SquareRatioTest: '+str(tagRatio) +' ID: '+str(self.__currentIdx)
                    if ratioSqualetoCircleRatio > 0.7:
                        '''
                        sanity check
                        '''
                        if (float(yMax) - float(yMin)) > 0:
                                '''
                                determine tag ratio
                                '''
                                tagRatio=(float(xMax) - float(xMin)) / (float(yMax) - float(yMin))
                                print 'Circle Ratio : '+str(tagRatio) +' ID: '+str(self.__currentIdx)
                                ratio.append(np.abs(1-(float(xMax) - float(xMin)) / (float(yMax) - float(yMin))))
                        else: ratio.append(1000)  
                    else: ratio.append(1000) 
                else: ratio.append(1000)  
            else: ratio.append(1000) 
        
    
        rect = np.min(ratio)
        rectIdx = list(ratio).index(rect)
        
        xMin = np.min(self.__compsX[rectIdx])
        xMax = np.max(self.__compsX[rectIdx])
        yMin = np.min(self.__compsY[rectIdx])
        yMax = np.max(self.__compsY[rectIdx])
       
        
        print 'Circle Ratio: '+str(rect)
        
        idx = np.where(labeled==rectIdx)
        '''
        bounding box
        '''
        iMin=np.min(idx[0])
        jMin=np.min(idx[1])
        iMax=np.max(idx[0])
        jMax=np.max(idx[1])
        
        sel = labeled != rectIdx
        labeled[sel]=0
        
        if rect > 0.2: 
            print 'Error: No circle detectable'
            rect=1
            rectIdx=0

        return rectIdx, rect, float(xMax) - float(xMin), float(yMax) - float(yMin),labeled[iMin:iMax, jMin:jMax]
                    

    def findRoot(self, labeledToCopy):
        print 'searching rootstock'
        labeled=labeledToCopy.copy()
        h,w=np.shape(labeled)
        found=False
        idx1=0
        count=0
        '''
        We keep this piece of debug code, because the problem occurs in 1 of 10,000 images. Perhaps we understand it one day.
        '''
        while found==False:
            idx1 = np.argmax(self.__labelHist)
            idx = np.where(labeled==idx1)
            if (np.max(idx[0])+1) == w and (np.max(idx[1])+1)==h and (np.min(idx[0])) == 0 and (np.min(idx[1]))==0:
                if count < len(self.__labelHist): 
                    found=False
                    count+=1
                else: found=True
                print 'Only 1 background component that is smaller than the foreground ??? Probably a bug in the Masking routine'
            else:
                found=True
    
        '''
        bounding box
        '''
        iMin=np.min(idx[0])
        jMin=np.min(idx[1])
        iMax=np.max(idx[0])
        jMax=np.max(idx[1])
       
        print 'xMin and xMax of Root Crown: '+str(iMin)+' '+str(iMax)
        print 'yMin and yMax of Root Crown: '+str(jMin)+' '+str(jMax)
        
        
        return idx1,idx,iMax,iMin,jMin,jMax
    
    def correctForStem(self,labeledToCopy,excludeIdx,left,right,bottom,top,rootIdx,rootIdxList):
        '''
        We loop through detected objects to identify them. During recognition the labeled image has to be made free of noise.
        Therefore we have to copy it.
        '''
        print 'checking for stem part'
        labeled=labeledToCopy.copy()
        idx2=-1
        counter=0
        again=True
        while again==True:
            counter+=1
            if counter>30: break
            print 'stem part loop'
            again=False
            nr_objects=len(self.__labelHist)
            if nr_objects>1:     
                comp = np.argmax(self.__labelHist)
                count=0
                if comp in excludeIdx:
                    self.__labelHist[comp] = 0      
                    comp = np.argmax(self.__labelHist)
                    if count>len(excludeIdx):
                        print 'Error: Image is not usable'
                        idx2=-1
                        break
                    else: count +=1
                idx2 = comp
                
            else:
                idx2=-1
            
            labeled=labeledToCopy.copy()
            idx = np.where(labeled==idx2)
            '''
            bounding box
            '''
            try:
                iMin=np.min(idx[0])
                jMin=np.min(idx[1])
                iMax=np.max(idx[0])
                jMax=np.max(idx[1])
                
                print 'xMin and xMax of stem part: '+str(iMin)+' '+str(iMax)
                print 'yMin and yMax of stem part: '+str(jMin)+' '+str(jMax)
                print 'yMax of crown: '+str(top)
                print 'xMin of crown: '+str(bottom)
                sel = labeled != idx2 
                labeled[sel]=0
                nonZ=len(idx[0])
                boundingBoxSize=(iMax-iMin)*(jMax-jMin)
                zeros=boundingBoxSize-nonZ
                ratio=float(zeros)/float(nonZ)
                print 'ratio: '+str(ratio)
                     
                if counter>=nr_objects:
                    again=False

            except:
                again=True
        nrOfObjPart=0
        if (right)>iMax*0.9: 
            
            sel = labeled == (idx2)
            self.__labelHist[idx2] = 0
            self.__labelHist[rootIdx]=0
            imgReturn=np.zeros_like(labeled)

            imgReturn[sel]=1
            imgReturn[rootIdxList]=1
            rep=int(np.fabs(iMax*0.9-right*1.1))
            for i in range(rep):
                imgReturn[iMax*0.9:right*1.1,jMin:jMax]=scipy.ndimage.binary_dilation(imgReturn[iMax*0.9:right*1.1,jMin:jMax]).astype(np.int)
                imgLabel,nrOfObjPart=scipy.ndimage.label(imgReturn[iMin:left,bottom:top])
                print 'nrOfObj = '+str(nrOfObjPart)
                if nrOfObjPart == 1:
                    break

            if nrOfObjPart ==1:
                return imgReturn[iMin:left,bottom:top]
            else:
                self.__labelHist[rootIdx]=0 
                imgReturn=np.zeros_like(labeled)
                imgReturn[rootIdxList]=1
                return imgReturn[right:left,bottom:top]
        else:
            self.__labelHist[rootIdx]=0 
            imgReturn=np.zeros_like(labeled)
            imgReturn[rootIdxList]=1
            return imgReturn[right:left,bottom:top]
    
    def findExcisedRoot(self, labeledToCopy,excludeIdx,minOfCrown,maxOfCrown):

        '''
        We loop through detected objects to identify them. During recognition the labeled image has to be made free of noise.
        Therefore we have to copy it.
        '''
        print 'searching excised root'
        labeled=labeledToCopy.copy()
        w,h=np.shape(labeled)
        idx2=-1
        counter=0
        again=True

        while again==True:
            counter+=1
            if counter>30: break
            print 'excised root loop'
            again=False
            nr_objects=len(self.__labelHist)
            if nr_objects>1:     
                comp = np.argmax(self.__labelHist)
                count=0
                if comp in excludeIdx:
                    
                    self.__labelHist[comp] = 0      
                    comp = np.argmax(self.__labelHist)
                    if count>len(excludeIdx):
                        print 'Error: Image is not usable'
                        idx2=-1
                        break
                    else: count +=1
                idx2 = comp
                self.__labelHist[idx2] = 0
            else:
                idx2=-1
            
            labeled=labeledToCopy.copy()
            idx = np.where(labeled==idx2)
            
            '''
            bounding box
            '''
            try:
                iMin=np.min(idx[0])
                jMin=np.min(idx[1])
                iMax=np.max(idx[0])
                jMax=np.max(idx[1])
                print 'xMin and xMax of Excised Root: '+str(jMin)+' '+str(jMax)
                print 'yMin and yMax of Excised Root: '+str(iMin)+' '+str(iMax)
                print 'xMax of crown: '+str(maxOfCrown)
                print 'xMin of crown: '+str(minOfCrown)
                sel = labeled != idx2
                labeled[sel]=0
                nonZ=len(idx[0])
                boundingBoxSize=(iMax-iMin)*(jMax-jMin)
                zeros=boundingBoxSize-nonZ
                ratio=float(zeros)/float(nonZ)
                print 'ratio: '+str(ratio)

                if counter>=nr_objects:
                    again=False

            except:
                again=True
        sel = labeled == idx2 
        labeled[sel]=255     
        return idx2,labeled[iMin:iMax, jMin:jMax],(iMax+iMin)/2,(jMax+jMin)/2,
        
    def findTag(self, labeled, imgBinary, img, rect_ratio=0.33):
        print 'searching tag'

        ratio = []
        for i in range(len(self.__compsX)):
            if self.__labelHist[i] > 0:
                xMin = np.min(self.__compsX[i])
                xMax = np.max(self.__compsX[i])
                yMin = np.min(self.__compsY[i])
                yMax = np.max(self.__compsY[i])
                '''
                The Tag should cover at least 0.5% of the picture
                '''
                if float((xMax-xMin)*(yMax-yMin))/float(self.__h*self.__w) >0.005:
                    '''
                    The Tag should have more length then height
                    '''
                    if (xMax-xMin)>=2*(yMax-yMin):
                        '''
                        The Tag should be in the upper half of the image
                        '''
                        if yMin< (self.__h*0.5): 
                            tagRatio=(float(xMax) - float(xMin)) / (float(yMax) - float(yMin))
                            print 'TagRatio detected: '+str(tagRatio)+' ID: '+str(self.__currentIdx)
                            ratio.append((float(xMax) - float(xMin)) / (float(yMax) - float(yMin)))
                        else: ratio.append(-1)
                    else: ratio.append(-1)
                else: ratio.append(-1)  
            else: ratio.append(-1) 
        
    
        rect = np.max(ratio)
        if rect >=0:
            rectIdx = list(ratio).index(rect)
            xMin = np.min(self.__compsX[rectIdx])
            xMax = np.max(self.__compsX[rectIdx])
            yMin = np.min(self.__compsY[rectIdx])
            yMax = np.max(self.__compsY[rectIdx])
        else: 
            xMin=xMax=yMin=yMax = 0
            rectIdx =-1

        print 'Tag Ratio: '+str(rect)
        if rect ==-1: 
            rectIdx=-1
            iMin=iMax=jMin=jMax=0
        else:
            idx = np.where(labeled==rectIdx)
            '''
            bounding box
            '''
            iMin=np.min(idx[0])
            jMin=np.min(idx[1])
            iMax=np.max(idx[0])
            jMax=np.max(idx[1])
        
        sel = labeled != rectIdx
        if rect>=0: 
            labeled[sel]=0
            try:
                print 'Check for text'
                tagText=ocr.getTextFromImage(img[iMin+self.__tagCrop:iMax-self.__tagCrop, jMin+self.__tagCrop:jMax-self.__tagCrop],self.__io.getHomePath(),str(self.__id))
            except:
                tagText='Tag text extraction Failed'
                pass
            '''
            This order prefers the barcode over the text reader
            '''
            try: 
                print 'Check for barcode'
                tagCode=ocr.getCodeFromImage(img[iMin+self.__tagCrop:iMax-self.__tagCrop, jMin+self.__tagCrop:jMax-self.__tagCrop],self.__io.getHomePath())
                if len(tagCode) >2:
                    tagText=tagCode[8:len(tagCode)-1]
                    print tagText
                else:
                    print 'bar code to short: '+tagCode
            except:
                pass
                
            self.__io.writeServerFile('dirt_out.csv',self.__io.getHomePath()+'/Mask/'+self.__io.getFileName()+'Tag.png,' +str(self.__io.getID())+',0')
            
            
        else:
            tagText= 'No label found'
        
        print 'return value rectIdx: '+str(rectIdx)
        return rectIdx, rect, float(xMax) - float(xMin), float(yMax) - float(yMin),imgBinary[iMin+self.__tagCrop:iMax-self.__tagCrop, jMin+self.__tagCrop:jMax-self.__tagCrop],tagText

    
