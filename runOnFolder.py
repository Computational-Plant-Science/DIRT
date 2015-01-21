'''
Created on Feb 12, 2014

@author: Alexander Bucksch
'''

'''
# external library imports
'''
import outputCrawler as oc
import time
'''
# python standard imports
'''
import os
import sys
import multiprocessing
import subprocess


def calculate(args):
    try:
        return subprocess.call(args)
    except:
        print "ERROR in File: "+str(args[2])

if __name__ == '__main__':
    print os.getcwd()
    startT=time.time()
    dir=sys.argv[1]  
    seg=sys.argv[2]  
    files=os.listdir(dir)
    pool = multiprocessing.Pool(processes=8)
    args=[]
    for idx,i in enumerate(files): 
        if i !='.DS_Store':
            if os.path.isfile(dir+i):
                if os.path.isdir(dir+str(idx))==False:
                    args.append(['python', os.getcwd()+'/main.py', dir+str(i),str(idx),seg, '1', '1', '1', '0.0', '0', '0', dir, './traits.csv'])
    r = pool.map_async(calculate, args)
    r.wait() # Wait on the results
    print 'All files done in '+str(time.time()-startT)+'s !'
    print 'Collecting results'
    oc.combineOutput(dir) 
    print 'Results written to '+dir+'outputAll.csv'
