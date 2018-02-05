'''
Created on May 21, 2013

@author: Alexander Bucksch
'''

'''
# standard python imports
'''
import os
import csv

def combineOutput(dirName):

    ''' save file'''
    
    ''' remove the old output '''
    try:
        os.remove(os.path.join(dirName, 'outputAll.csv'))
    except:
        pass
    try:
        os.remove(os.path.join(dirName, 'crawler.out'))
    except:
        pass
    files = os.listdir(dirName)
    directories=[]
    for i in files:
        if os.path.isdir(os.path.join(dirName, i)):
            directories.append(i)
            
    ''' copy header'''
    print directories 
    with open(os.path.join(dirName, directories[0], 'output.csv'), 'U') as csvfile:
        filedata = csv.reader(csvfile)
        rows = filedata.next()
    
        with open(os.path.join(dirName, 'outputAll.csv'), 'w') as f:
            filewriter = csv.writer(f)
            filewriter.writerow(rows)
    ''' append data'''
    countOK = 0
    countBAD = 0
    badFolder = []
    for f in directories:
        try:
            with open(os.path.join(dirName , f , 'output.csv'), 'U') as csvfile:
                countOK += 1
                filedata = csv.reader(csvfile)
                with open(os.path.join(dirName, 'outputAll.csv'), 'a+') as of:
                    filewriter=csv.writer(of)
                    rows=filedata.next()
                    rows=filedata.next()
                    filewriter.writerow(rows)
        except:
            countBAD += 1
            badFolder.append(f)
            # Open a file
    with open(os.path.join(dirName, "crawler.out"), "wb") as fo:
        fo.write( str(countOK) + ' images are processed and ' +str(countBAD) + ' images failed: \n' + str(badFolder))
        
    print str(countOK)+' images are processed and '+str(countBAD)+' images failed: \n'+str(badFolder)
            