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

    def __init__(self, homePath=None, name=None, ID=None, plots=True):
        '''
        Constructor
        '''
        self.__plots = plots
        self.__id = ID
        self.__currentID = None
        self.__path = homePath
        self.__name = name
        self.__serverPath = None
        self.__parameters = ['Image ID', 'Image name', 'Failed',
                             'Experiment number',
                             'circle ratio',
                             'x pixel',
                             'y pixel',
                             'xScale',
                             'yScale', 'computation time', 'Skeleton Vertices']

    def setServerPath(self, p):
        self.__serverPath = p

    def setFileName(self, name):
        self.__name = name

    def getFileName(self):
        return self.__name

    def getHomePath(self):
        return self.__path

    def getID(self):
        return self.__id

    def getCurrentID(self):
        return self.__currentID

    def setidIdx(self, idx):
        self.__currentID = idx

    def setHomePath(self, homePath):
        self.__path = homePath

    def scanDir(self, directory=0):

        scanPath = ''
        if directory == 0:
            scanPath = self.__path
        else:
            scanPath = directory
        files = []
        print(os.getcwd())
        listing = os.listdir(scanPath)
        for infile in listing:
            if os.path.isdir(scanPath + infile) == True:
                print(infile + ' is not a file')
            elif infile[0] == '.':
                print(infile + ' is not a file')
            elif infile[len(infile) - 4:] == '.ini':
                print(infile + ' is not a file')
            elif infile[len(infile) - 4:] == '.csv':
                print(infile + ' is not a file')
            else:
                print("current file is: " + infile)
                files.append(scanPath + '/' + infile)
        return files

    def writeServerFile(self, serverFile, string):
        path = self.__serverPath
        print("server file (working directory): " + str(os.getcwd()))
        print("server file (relative): " + str(path))
        try:
            if os.path.isfile(path + serverFile):
                fout = open(path + serverFile, "a")
            else:
                raise
        except:
            fout = open(path + serverFile, "w")
            fout.write('# File path, image id, type')
            fout.write('\n')

        fout.write(string)
        fout.write('\n')
        fout.close()

    def writeRunFile(self, runfile, string):

        path = self.__serverPath + '/'
        try:
            if os.path.isfile(path + 'dirtRun.csv'):
                fout = open(path + 'dirtRun.csv', "a")
            else:
                raise
        except:
            fout = open(path + 'dirtRun.csv', "w")
            fout.write('# File path, image id')
            fout.write('\n')

        fout.write(runfile + ', ' + string)
        fout.write('\n')
        fout.close()

    def writeFile(self, para, traitsCrown, traitDict, all=False):
        print("output directory: " + self.__path + "/output.csv")
        try:
            if os.path.isfile(self.__path + "/output.csv"):
                fout = open(self.__path + "/output.csv", "a")
            else:
                raise
        except:
            fout = open(self.__path + "/output.csv", "w")
            for i in self.__parameters:
                fout.write(str(i) + ',')
            for k, v in traitDict.items():
                if all == False:
                    if v == True:
                        if k in traitsCrown:
                            fout.write(str(k) + ',')
                else:
                    fout.write(str(k) + ',')

            fout.write('\n')

        for i in para:
            fout.write(str(i) + ',')
        for k, v in traitDict.items():
            if v == True and k in traitsCrown:
                fout.write(str(traitsCrown[k]) + ',')
            elif all == True:
                fout.write(str(' ,'))

        fout.write('\n')
        fout.close()

    def saveArray(self, arr, name):
        if self.__plots == False: return 0
        try:
            np.savetxt(name + '.gz', arr, delimiter=',')
            try:
                self.writeServerFile('dirt_out.csv', os.getcwd() + name[1:] + '.gz' + ',' + str(self.__id) + ',1')
            except:
                raise
        except:
            raise
