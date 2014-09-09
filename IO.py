'''
IO.py

We handle the input and output of the software in this module

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

'''
# standard python imports
'''
import os


class IO(object):
    '''
    classdocs
    '''
    
    __instance = None

    ## Class used with this Python singleton design pattern
    #  @todo Add all variables, and methods needed for the Singleton class below
    class Singleton:
        def __init__(self):
            ## a foo class variable
            self.__path = None
            self.__name = None
    
    def __init__(self,homePath,name=None,ID=None):
        '''
        Constructor
        '''
        self.__id=ID
        self.__currentID=None
        self.__path = homePath
        self.__name = name
        self.__serverPath=None
        self.__traitsCrown=['stemDia','avg. Root Density', 
                 'STA range', 
                 'STA dom. angle 1', 
                 'STA dom. angle 2',' STA 25% 1','STA 25% 2','STA 50% 1','STA 50% 2','STA 75% 1','STA 75% 2','STA 90% 1','STA 90% 2',
                 'RTA dom. angle 1', 
                 'RTA dom. angle 2',
                 'STA min.',
                 'STA max.',
                 'STA median',
                 'RTA range',
                 'RTA min.',
                 'RTA max.',
                 'RTA median',
                 'Nr. of RTPs', 
                 'TD median',
                 'TD mean',
                 'DD90 max.','median Width','max. Width','D10','D20','D30','D40','D50','D60','D70','D80','D90','DS10','DS20','DS30','DS40','DS50','DS60','DS70','DS80','DS90',
                 'spatial root distr. X','spatial root ditr. Y','CPD 25','CPD 50','CPD 75','CPD 90']

        self.__traitsLateral=['avg. lateral length',
                       'nodal root path length',
                       'lateral branching frequency',
                       'avg. nodal root diameter',
                       'nodal root shape',
                       'lateral avg. angle',
                       'lateral angular range',
                       'lateral min. angle',
                       'lateral max. angle',
                       'dist to first lateral',
                       'median lateral diameter','mean lateral diameter'
                       ]
        
        self.__parameters=['Image ID','Image name','Failed',
                           'Experiment number',
                    'circle ratio',
                    'x pixel',
                    'y pixel',
                    'xScale',
                    'yScale','computation time','Skeleton Vertices']

            
    def setServerPath(self,p):
        self.__serverPath=p    
    def setFileName(self,name):
        self.__name = name
    def getFileName(self):
        return self.__name    
    def getHomePath(self):
        return self.__path
    def getID(self):
        return self.__id
    def getCurrentID(self):
        return self.__currentID
    def setidIdx(self,idx):
        self.__currentID=idx
    def setHomePath(self,homePath): 
        self.__path = homePath
    
    def scanDir(self,directory=0):

        scanPath=''
        if directory==0: scanPath=self.__path
        else: scanPath = directory
        files = []
        print os.getcwd()
        listing = os.listdir(scanPath)
        for infile in listing:
            if os.path.isdir(scanPath + infile) == True:
                print infile + ' is not a file'
            elif infile[0]=='.':
                print infile + ' is not a file'
            elif infile[len(infile)-4:]=='.ini':
                print infile + ' is not a file'
            elif infile[len(infile)-4:]=='.csv':
                print infile + ' is not a file'
            else:
                print "current file is: " + infile
                files.append(scanPath+'/'+infile) 
        return files
    
    def writeServerFile(self,serverFile,string):
        pass
    
    def writeFile(self,para,traitsCrown,traitsLateral):
        try:
            if os.path.isfile(self.__serverPath+"/output.csv"):
                fout=open(self.__serverPath+"/output.csv", "a")
            else: raise
        except:
            fout = open(self.__serverPath+"/output.csv", "w")
            for i in self.__parameters:
                    fout.write(str(i)+',')
            for i in self.__traitsCrown:
                    fout.write(str(i)+',')
            for i in self.__traitsLateral:
                    fout.write(str(i)+',')
            fout.write('\n')
            

    
        for i in para:
                    fout.write(str(i)+',')
        for i in traitsCrown:
                    fout.write(str(i)+',')
        for i in traitsLateral:
                    fout.write(str(i)+',')
        fout.write('\n')
        fout.close()
        
    def saveArray(self,arr,name):
        pass
        
        
