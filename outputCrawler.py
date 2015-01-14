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
        os.remove(dirName+'outputAll.csv')
    except:
        pass
    try:
        os.remove(dirName+'crawler.out')
    except:
        pass
    files=os.listdir(dirName)
    directories=[]
    for i in files:
        if os.path.isdir(dirName+i)==True:
            directories.append(i)
            
    ''' copy header'''
    with open (dirName+directories[0]+'/output.csv','U') as csvfile:
        filedata= csv.reader(csvfile)
        rows=filedata.next()
    
        with open(dirName+'outputAll.csv', 'w') as f:
            filewriter=csv.writer(f)
            filewriter.writerow(rows)
    ''' append data'''
    countOK=0
    countBAD=0
    badFolder=[]
    for f in directories:
        try:
            with open (dirName+f+'/output.csv','U') as csvfile:
                countOK+=1
                filedata= csv.reader(csvfile)
                with open(dirName+'outputAll.csv', 'a+') as f:
                    filewriter=csv.writer(f)
                    rows=filedata.next()
                    rows=filedata.next()
                    filewriter.writerow(rows)
        except:
            countBAD+=1
            badFolder.append(f)
            # Open a file
    fo = open(dirName+"crawler.out", "wb")
    fo.write( str(countOK)+' images are processed and '+str(countBAD)+' images failed: \n'+str(badFolder))
    # Close opend file
    fo.close()
    print str(countOK)+' images are processed and '+str(countBAD)+' images failed: \n'+str(badFolder)
            
if __name__ == '__main__':
    combineOutput('/Users/koalaspirit/Documents/DIRTTestset/')