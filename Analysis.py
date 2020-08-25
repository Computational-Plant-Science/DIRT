'''
Analysis.py

The analysis module for DIRT. Most of the traits are computed here.

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
# internal library imports
'''
import kmeans as km
import ransac
'''
# external library imports
'''
import numpy as np
import warnings
warnings.simplefilter('ignore', np.RankWarning) #suppress the warning
from scipy import polyfit, polyval
import scipy.stats
import scipy.misc
import scipy.interpolate
from graph_tool.all import *
import graph_tool.topology as gt

class Analysis(object):
    '''
    classdocs
    '''


    def __init__(self,io,scale):
        '''
        Constructor
        '''
        self.__io = io
        self.__id=io.getID()
        self.__currentIdx=io.getCurrentID()
        self.__scale=scale
    

    def findHistoPeaks(self,ang):
        try:
            pdf, _ =np.histogram(ang,bins=9, range=(0, 90))
            maximum=[]
            for i in range(0,len(pdf)-1):
                if pdf[i-1] <= pdf[i] or i==0:
                    if pdf[i+1] <= pdf[i] or i==len(pdf)-1:
                        maximum.append(i)
            
            pdfSort=np.argsort(pdf)
            strongestMax = []
            for i in pdfSort[::-1]:
                if i in maximum:
                    strongestMax.append(i)
                if len(strongestMax)==2:
                    break
            avgAng=[]
            avgResult=[]
            for i in strongestMax:
                avgAng=[]
                for j in ang:
                    if j >= i*10:
                        if j <= (i+1)*10:
                            avgAng.append(j)
                avgResult.append(np.mean(avgAng))
            if (len(avgResult)>1):
                return np.max(avgResult),np.min(avgResult)
            else:
                return avgResult[0],-1
            
        except:
            print 'WARNING: no histo peaks found: Analysis.findHistoPeaks'
            print pdf
            return -1,-1

    def smooth(self,x,window_len=11,window='hanning'):
        """smooth the data using a window with requested size.
        
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
            raise ValueError, "smooth only accepts 1 dimension arrays."
    
        if x.size < window_len:
            raise ValueError, "Input vector needs to be bigger than window size."
    
    
        if window_len<3:
            return x
    
    
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    
    
        s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
        #print(len(s))
        if window == 'flat': #moving average
            w=np.ones(window_len,'d')
        else:
            w=eval('np.'+window+'(window_len)')
    
        y=np.convolve(w/w.sum(),s,mode='valid')
        return y

    def filterPathDiameters(self,path,G,scale):
        #remove diameters around branching points in range of the branching point diameter
        vprop=G.vertex_properties["vp"]
        for i in range(len(path)):
            count = 0
            for _ in path[i].out_neighbours():
                count+=1
                if count >2:
                    break
            if count>2:
                for j in range(int(vprop[path[i]]['diameter']/self.__scale)):
                    if j >20: break
                    if i-j >0: vprop[path[i-j]]['diameter']=0
                    if i+j <len(path):vprop[path[i+j]]['diameter']=0
                
                vprop[path[i]]['diameter']=0
        return G

    def filterPath(self,path,G):
        #remove diameters around branching points in range of the branching point diameter
        delDia=G.node[path[0]]['diameter']
        
        return path[3::]
                    
    def getLateralLength(self,pathList,thickestPath,G,counter=None):
        vprop=G.vertex_properties["vp"]
        lengthArr=[]
        x=[]
        if len(pathList)>0:
            for i in pathList:
                length = len(i)
                if len(i)>0:
                    lengthArr.append(length)  
                    x.append(vprop[G.vertex(i[0])]['imgIdx'][0])  
                else:   
                    lengthArr.append(-1)  
                    x.append(-1)  
        avgLength=np.average(lengthArr)
        
        print 'avg. Length: ' + str(avgLength)
        try:
            self.__io.saveArray(lengthArr,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_LengthHist')

            f2 = scipy.interpolate.interp1d(x, lengthArr, kind='cubic')
            self.__io.saveArray(x,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_LengthX')
            self.__io.saveArray(f2(x),self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_LengthY')
            
        except:
            pass
        print "Avg.Length Scale Nodalroot"+str(self.__scale)
        return avgLength*self.__scale
        
    def getLateralLengthRTP(self,RTP,img,counter=None):
        lengthArr=[]
        imgDebug=img.astype(np.uint8)
        for i in RTP:
            length = len(i)
            for j in i:
                imgDebug[j[1]][j[0]]=255
                lengthArr.append(length)
            scipy.misc.imsave(self.__io.getHomePath()+'Result/' +str(j[1])+' '+str(j[0])+ 'DebugThickP.png', imgDebug)  
            self.__io.writeServerFile(self.__io.getHomePath(), 'dirt_out.csv',self.__io.getHomePath()+'Result/'+self.__io.getFileName()+'DebugThickP.png,' +str(self.__id[self.__currentIdx])+',0')   
        avgLength=np.average(lengthArr)
        print 'avg. Length: ' + str(avgLength)
        
        self.__io.saveArray(lengthArr,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_LateralLengthHisto')
        self.__io.saveArray(range(0,len(lengthArr)),self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_LateralLengthX')
        self.__io.saveArray(lengthArr,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_LateralLengthY')
        
    def getSymmetry(self,rtps,G):
        vprop=G.vertex_properties["vp"]
        #calculate bounding box
        xMax=0
        xMin=5000000000
        yMax=0
        yMin=5000000000
        sumX=0
        sumY=0
        count=0.0
        for r in rtps:
            for i in r:
                v=G.vertex(i) 
                if vprop[v]['coord'][0] > xMax:
                    xMax=vprop[v]['coord'][0]
                if vprop[v]['coord'][0] < xMin:
                    xMin=vprop[v]['coord'][0]
                if vprop[v]['coord'][1] > yMax:
                    yMax=vprop[v]['coord'][1]
                if vprop[v]['coord'][1] < yMin:
                    yMin=vprop[v]['coord'][1]
                sumX+=vprop[v]['coord'][0]
                sumY+=vprop[v]['coord'][1]
                count+=1.0
        if count>0:
            avgX=sumX/count
            avgY=sumY/count
            avgBoxX=(xMax-xMin)/2
            avgBoxY=(yMin-yMax)/2
            vecX=avgBoxX-avgX
            vecY=avgBoxY-avgY
            vecSym=[vecX,vecY]
        else:
            vecSym=[np.nan,np.nan]
            
        return vecSym
    
    def getLateralAngles(self,thickestPath,lat,corrBranchpts,G,counter=None):
        # calculates the angle between the hypocotyle and all rtps
        angles=[]
        tangents=self.filterRTPTangent(thickestPath,lat,corrBranchpts)
        for i in range(len(lat)):
            ang=self.getAngleBetweenPaths(tangents[i],lat[i],G)
            angles.append(ang)
        try:
            self.__io.saveArray(angles,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_AngleHisto')
            minAngle=np.min(angles)
            maxAngle=np.max(angles)
            angRange=maxAngle-minAngle
            avgAngle=np.average(angles)

        except: 
            minAngle=-1
            maxAngle=-1
            angRange=-1
            avgAngle=-1
        return angRange,avgAngle,minAngle,maxAngle,angles
                
    def getBranchingfrequencyAlongSinglePath(self,rtps,path):
        bp=[]
        for i in rtps:
            if i[0]!=path[0]:
                bp.append(i[0])
        
        bpUnique=np.unique(bp)    
        try:
            branchFreqency=float(len(path))/float(len(bpUnique))
        except:
            branchFreqency=-1
        print 'Branching Frequency in given unit: ' +str(branchFreqency*self.__scale)
        return branchFreqency*self.__scale
        
    def getDiametersAlongSinglePath(self,path,G,scale,counter=None):
        vprop=G.vertex_properties["vp"]

        x=[]
        y=[]

        length=0
        for i in path:
            length+=1
            if vprop[i]['diameter'] > 0:
                x.append(length)
                y.append(vprop[i]['diameter'])
        coeffs=polyfit(x,y,1)

        besty =  polyval ( coeffs ,    x)
        
        self.__io.saveArray(x,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_DiameterX')
        self.__io.saveArray(y,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_DiameterY')

        avgDiameter=np.average(y)
        self.__io.saveArray(y,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_DiameterHisto')
        return avgDiameter,coeffs[0]
        
    def getDiameterQuantilesAlongSinglePath(self,path,G,counter=None):
        
        G=self.filterPathDiameters(path, G,self.__scale)
        x=[]
        y=[]
        length=0
        vprop=G.vertex_properties["vp"]
        for i in path:
            length+=1
            if vprop[i]['diameter'] > 0:
                x.append(length)
                y.append(vprop[i]['diameter'])
        coeffs=polyfit(x,y,1)

        besty =  polyval ( coeffs ,    x)
        
        self.__io.saveArray(x,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_DiameterX')
        self.__io.saveArray(y,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_DiameterY')

        l=len(y)-1
        l25=int(l*0.25)
        l50=int(l*0.5)
        l75=int(l*0.75)
        l90=int(l*0.90)
        
        d25=np.average(y[:l25])
        d50=np.average(y[l25:l50])
        d75=np.average(y[l50:l75])
        d90=np.average(y[l90:])

        
        self.__io.saveArray(y,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_DiameterHistoTP')
        
        return d25,d50,d75,d90
            
    def getWidthOverHeight(self,img2,xScale,yScale):
        # We compute here all mask based traits at once (long lives Spagetti code :-) )
        if len(img2) > 0: 
            print 'IMG OK'
        else:
            print 'IMAGE IS NOT ACCESSIBLE'
            return 'nan','nan','nan','nan','nan','nan','nan','nan',['nan']*9,['nan']*9,'nan'
        h, w = np.shape(img2)
        xx=[]
        yy=[]
        blackArr=[]
        sizeCount=0
        densityArray=[]
        #compute density value
        for i in range(h):
            white=0
            black=1
            start = 0
            end = 0
            x = i / w
            
            idx = np.where(img2[i]>0)
            try:
                start=idx[0][0]
                white=len(idx[0])
                end=idx[0][len(idx[0])-1] 
                black=(end-start)-white
                sizeCount+=white
                width=float(end-start)
                xx.append(float(i)*yScale)
                if black==0: black=1
                blackArr.append(black)
                yy.append(float(width)*xScale)
                
                normWhite=1
                normBlack=1

                if w > 0:
                    normWhite=float(white)/float(w)
                    if normWhite==1.: normWhite=0.
                    if black >0: normBlack=float(black)/float(w)
                    else: normBlack=1.

                    if normWhite>0.: densityArray.append(float(normWhite)/float(normBlack))
            except:
                densityArray.append(float(0.0))
                print 'empty image line in crown file -> placed 0. as density for this line'
                pass
            
        rootDensity=np.average(densityArray)
        print 'Avg. Root density: ' + str(rootDensity)
        
        ysmooth = yy
        smoothRegion=15
        xxNorm=np.array(xx)/np.max(xx)
        #changed to 1% for experiment 
        tenPercent=float(len(yy))*0.01
        # retrieve stem diameter as the average if the distance field in the first 10%
        try:
         stemDia=np.median(xx[20:int(tenPercent)])
        except:
         stemDia=-1
        # compute a simple angle at top and bottom along the outline for monocots (note this is more noisy than the D10 or D20 values that are more robust)
        try:
         #changed to 1 to 3%
         angleSimple=self.getAngleToXAxXY(xx[int(tenPercent):int(tenPercent)*3], yy[int(tenPercent):int(tenPercent)*3], ransacFitting=True)

         #changed to 3-5%
         angleSimpleBottom=self.getAngleToXAxXY(xx[int(tenPercent)*3:int(tenPercent)*5], yy[int(tenPercent)*3:int(tenPercent)*5], ransacFitting=True)
        except:
         angleSimple=-1
         angleSimpleBottom=-1
        print 'Stem Diameter: '+str(stemDia)
        print 'Root Top Angle: '+str(angleSimple)
        print 'Root Bottom Angle: '+str(angleSimpleBottom)
        
        #smooth the noisy data
        for i in range(smoothRegion,len(yy)-smoothRegion):
            tmp=[]
            for j in range(smoothRegion):
                tmp.append(yy[i+j])
                tmp.append(yy[i-j])
            ysmooth[i]=np.median(tmp)   
        self.__io.saveArray(xxNorm,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_HeightWidthX') 
        self.__io.saveArray(ysmooth,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_HeightWidthY') 
        
        # compute the width parameters of the root
        medianWidth=np.median(ysmooth)
        maxWidth=np.max(ysmooth)
        print 'Median. Root Width: ' + str(medianWidth)
        print 'Max. Root Width: ' + str(maxWidth)
        #compute the cummulative width profile
        ysmoothCS=np.array(ysmooth).cumsum()
        ysmoothCS=ysmoothCS/np.max(ysmoothCS)
        
        #find index at x% width accumulation
        D=[]
        
        dD=0.1
        count=0.0
        for i in ysmoothCS:
            count+=1.0
            if i>dD:
                D.append(count)
                dD+=0.1
            if dD==1.0:
                break
        while len(D)<9:
            D.append(-1)
        
        
        # for drawing we need all DS valuse ofver the curve
        #Dslope=self.filterCDFTAngentSlope(xxNorm,ysmoothCS, range(len(ysmoothCS)),200)
        
        #save the data for plotting
        self.__io.saveArray(xxNorm,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_HeightWidthCSX') 
        self.__io.saveArray(ysmoothCS,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_HeightWidthCSY')
        #self.__io.saveArray(Dslope,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_HeightWidthDSX')
        #self.__io.saveArray(ysmoothCS,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_HeightWidthDSY')

        # for the final output we just need the 10 reveant slopes
        # Note, we adjusted the window size to 100. In the paper we used 20
        print '****** filterCDFTAngentSlope Input **********'
        print xxNorm
        print ysmoothCS
        print D[:len(D)-1]
        
        Dslope=self.filterCDFTAngentSlope(xxNorm,ysmoothCS, D[:len(D)-1],100)
        #convert D array counts to percentages saved in the output file
        D=np.array(D,dtype=float)/float(len(ysmoothCS))
        return rootDensity,medianWidth,maxWidth,D[0],D[1],D[2],D[3],D[4],D[5],D[6],D[7],D[8],Dslope[0],Dslope[1],Dslope[2],Dslope[3],Dslope[4],Dslope[5],Dslope[6],Dslope[7],Dslope[8],sizeCount*(xScale*yScale),stemDia,angleSimple,angleSimpleBottom

    def getLengthOfPath(self, path):
        length=len(path)*self.__scale
        return length
    
    def calculateAngleAtDist(self,thickestPath,lat,corrBranchpts,scale,G,atDist=20,counter=None):
        # calculates the angle between the hypocotyle and all rtps
        
        angelsAtDist=[]

        pxDist=int(atDist/scale)
        for i in range(len(lat)):
            if len(lat[i])>pxDist and G.vertex(corrBranchpts[i]) in thickestPath:
                angAtDist=self.getAngleToXAx(G,lat[i][:pxDist])
                angelsAtDist.append(angAtDist)
            
        self.__io.saveArray(angelsAtDist,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_AngleHistoAtDist')
        return np.mean(angelsAtDist)
    
    def calculateAngleQuantiles(self,thickestPath,lat,corrBranchpts,G,counter=None):
        
        angles25=[]
        angles50=[]
        angles75=[]
        angles90=[]
        
        for i in range(len(lat)):
            l=len(lat[i])
            l25=int(l*0.25)
            l50=int(l*0.5)
            l75=int(l*0.75)
            l90=int(l*0.90)
            
            ang25=self.getAngleToXAx(G,lat[i][:l25])
            ang50=self.getAngleToXAx(G,lat[i][:l50])
            ang75=self.getAngleToXAx(G,lat[i][:l75])
            ang90=self.getAngleToXAx(G,lat[i][:l90])
            angles25.append(ang25)
            angles50.append(ang50)
            angles75.append(ang75)
            angles90.append(ang90)
            
        self.__io.saveArray(angles25,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_AngleHisto25')
        self.__io.saveArray(angles50,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_AngleHisto50')
        self.__io.saveArray(angles75,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_AngleHisto75')
        self.__io.saveArray(angles90,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_AngleHisto90')
        
        return angles25,angles50,angles75,angles90

    def calculateAngles(self,thickestPath,lat,corrBranchpts,G,counter=None):
        # calculates the angle between the hypocotyle and all rtps
        angles=[]
        #tangents=self.filterRTPTangent(thickestPath,lat,corrBranchpts,window=45)
        #RTP=self.filterRTP(RTP,thickestPath,skel)
        for i in range(len(lat)):
            ang=self.getAngleToXAx(G,lat[i])
            angles.append(ang)
        self.__io.saveArray(angles,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_AngleHistoN')
#        f=p.figure()
#        ax = f.add_subplot(111)
#        p.hist(angles,bins=9, range=(0, 90))
#        p.title('Histogram of angles between hypocotyle and roots')
#        p.xlabel('angles in degree')
#        p.ylabel('# of angles')
        minAngle=np.min(angles)
        maxAngle=np.max(angles)
        angRange=maxAngle-minAngle
        avgAngle=np.median(angles)
#        p.text(0.3, 0.9,'angular range: '+str(angRange), ha='center', va='center', transform=ax.transAxes)
#        p.text(0.3, 0.85,' median angle: '+str(avgAngle), ha='center', va='center', transform=ax.transAxes)
#        p.savefig(self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_AngleHisto.png')
#        self.__io.writeServerFile('./', 'dirt_out.csv',self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_AngleHisto.png,' +str(self.__id[self.__currentIdx])+',1')
#        p.clf()
        return avgAngle,minAngle,maxAngle,angRange,angles

    def filterRTP(self,RTP,thickestpath,skel):
        RTPnew =[]
        for i in range(len(RTP)):
            RTPnew.append(frozenset(RTP[i]).difference(frozenset(thickestpath)))

            diameter=skel.node[i[0]]['diameter']

            for j in range(int(diameter)):
                i.pop(j)
                
        return RTPnew
    
    def filterCDFTAngentSlope(self,x,CDF,idx, window=20):
        #We estimate the slope at a point over a small region along the CDF
        tangents=[]
        for i in idx:
            tmpTangentX=[]
            tmpTangentY=[]
            for j in range(window):
                print 'appendices'
                print x[int(i)-j],CDF[int(i)-j]
                if j<i:
                    try:
                        tmpTangentX.append(x[int(i)-j])
                        tmpTangentY.append(CDF[int(i)-j])
                    except:
                        pass
                if j<len(CDF)-j:
                    try:
                        tmpTangentX.append(x[int(i)+j])
                        tmpTangentY.append(CDF[int(i)+j])
                    except: pass
            a,_=self.fitLineXY(tmpTangentX,tmpTangentY)
            tangents.append(a)
        return tangents
    
    def filterRTPTangent(self,thickestPath,lat,corrBranchpts, window=5):
        tangents=[]
        for i in corrBranchpts:
            tmpTangent=[]
            try:
                idx=thickestPath.index(i)
                
                for j in range(window):
                    if j<idx:
                        try:tmpTangent.append(thickestPath[idx-j])
                        except: pass
                    if j<len(thickestPath)-j:
                        try: tmpTangent.append(thickestPath[idx+j])
                        except: pass
            except: pass
            tangents.append(tmpTangent)
        return tangents
    
    def getAngleBetweenPaths(self, path1,path2,G,counter=None):
        #Helper function that calculates the angle between two given paths
        m1,_=self.fitLine(path1,G)
        m2,_=self.fitLine(path2,G)
        
        up=m1-m2
        low=1+(m1*m2)
        tanAlpha = np.fabs(up/low)
        alpha= (np.arctan(tanAlpha)*180)/np.pi
        return np.fabs(alpha)
    
    def getAngleToXAx(self, G,path2,counter=None):
        m2,_=self.fitLine(path2,G)
        alpha= (np.arctan(m2)*180)/np.pi
        return np.fabs(alpha)

    def getAngleToXAxXY(self, X,Y,counter=None,ransacFitting=False):
        m2,_=self.fitLineXY(X,Y,ransacFitting=ransacFitting)
        alpha= (np.arctan(m2)*180)/np.pi
        return np.fabs(alpha)
    
    def plotDiaRadius(self,paths,dia,thickestPath,nrOfClusters):
        print 'do the kmeans :-)'
        pts=[]
        for i in range(len(paths)):
            pts.append([paths[i],dia[i]])
        
        cl=km.kMeans(pts)
        c=cl.kmeans(nrOfClusters, 0.01)
        
        cX=[]
        cY=[]
        
        for nc in range(nrOfClusters):
            cX.append([])
            cY.append([])
        for nc in range(nrOfClusters):
            for i in c[nc]:
                cX[nc].append(i[0])
                cY[nc].append(i[1])
        for nc in range(nrOfClusters):
            self.__io.saveArray(cX[nc],self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_PathsDiaAx_'+str(nrOfClusters)+'_'+str(nc))
            self.__io.saveArray(cY[nc],self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_PathsDiaAy_'+str(nrOfClusters)+'_'+str(nc))

        if nrOfClusters ==2:
            if cY[0][0]>cY[1][0]:
                return cX[0],cY[0],cX[1],cY[1]
            else:
                return cX[1],cY[1],cX[0],cY[0]
        if nrOfClusters ==3:
            if cY[0][0]>cY[1][0] and cY[1][0]>cY[2][0]:
                return cX[0],cY[0],cX[1],cY[1],cX[2],cY[2]
            if cY[0][0]>cY[1][0] and cY[1][0]<cY[2][0]:
                return cX[0],cY[0],cX[2],cY[2],cX[1],cY[1]
            if cY[0][0]<cY[1][0] and cY[1][0]<cY[2][0]:
                return cX[2],cY[2],cX[1],cY[1],cX[0],cY[0]
            else:
                return cX[2],cY[2],cX[0],cY[0],cX[1],cY[1]
         
    def fitLine(self,path,G):
        #Simple line fit.
        vprop=G.vertex_properties["vp"]
        if len(path) == 0:
            return -1,-1
        X =[]
        Y= []
        for m in path:
            X.append(vprop[G.vertex(m)]['coord'][0])
            Y.append(vprop[G.vertex(m)]['coord'][1])
        (ar,br)=polyfit(X,Y,1)
        return ar,br
    def fitLineXY(self,X,Y, ransacFitting=False):
        #Simple line fit.
        if ransacFitting:
            try: 
                X,Y=ransac.ransacFit(X,Y)
            except: 
                print "ransac fitting failed. Using simple linear fitting"
        try:
		print '**** Polyfit Input *****'
		print X
		print Y
		(ar,br)=polyfit(X,Y,1)
	except:
		print "Fitting Failed (Analysis.py line 709)"
		ar=1.
		br=1.
        return ar,br

    def countRootsPerSegment(self, cAdv,cBas, cDiaAdv, cDiaBas):
        nrOfAdvRoots=0
        nrOfBasRoots=0
        firstSegmentRoots=cAdv[0]-cAdv[len(cAdv)-1]
        secondSegmentRoots=cBas[0]-cBas[len(cBas)-1]
        
        for idx,i in enumerate(cAdv):
            if cAdv[idx-1]-i > 0: #just a little treshhold to get rid of noise
                nrOfAdvRoots+=1
        for idx,i in enumerate(cBas):
            if cBas[idx-1]-i > 0: #just a little treshhold to get rid of noise
                nrOfBasRoots+=1
                
        return nrOfAdvRoots,nrOfBasRoots,firstSegmentRoots,secondSegmentRoots, np.mean(cDiaAdv),np.mean(cDiaBas)
            
    def RTPsOverDepth(self, centralPath, rtpSkel):
        
        vprop=rtpSkel.vertex_properties['vp']
        depth=[len(centralPath)]
        nrOfP=[len(centralPath)]
        fiftyPercentRtp=vprop[centralPath[0]]['nrOfPaths']/2
        fiftyPercentDrop=0
        for i in centralPath:
            nrOfP.append(vprop[i]['nrOfPaths'])
            depth.append(vprop[i]['coord'][1])
            if vprop[i]['nrOfPaths']<=fiftyPercentRtp:
                fiftyPercentDrop=vprop[i]['coord'][1]
                
        
        self.__io.saveArray(nrOfP,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_RTPDepthX')
        self.__io.saveArray(depth,self.__io.getHomePath()+'Plots/'+self.__io.getFileName()+'_RTPDepthY')
        
        return fiftyPercentDrop
    def anglesPerClusterAtDist(self, cAdv, cBas, rtpSkel, path, lat,corrBranchpts, scale, dist=20):
        # estimating the angles per cluster leads to distinguishing adventious roots from basal roots.
        
        vprop=rtpSkel.vertex_properties["vp"]
        minPathsAdv=np.min(cAdv)
        for idx,i in enumerate(path):
            if minPathsAdv>vprop[i]['nrOfPaths']:
                meanAdv=self.calculateAngleAtDist(path[:idx],lat,corrBranchpts,scale,rtpSkel,dist,None)
                meanBas=self.calculateAngleAtDist(path[idx:],lat,corrBranchpts,scale,rtpSkel,dist,None)
                break
        print 'Adv and Bas Angle at 2cm:'+str(meanAdv)+' , '+str(meanBas)        
        return meanAdv,meanBas
    
