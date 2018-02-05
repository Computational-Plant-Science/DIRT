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
    main_py_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'main.py'))
    startT=time.time()
    dir = os.path.abspath(sys.argv[1])
    work_dir = os.path.abspath(sys.argv[2])
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)        
    seg=sys.argv[3]  
    files=os.listdir(dir)
    pool = multiprocessing.Pool(processes=8)
    args=[]
    for idx,i in enumerate(files): 
        if i != '.DS_Store' and os.path.isfile(os.path.join(dir, i)):
            if not os.path.isdir(os.path.join(dir, str(idx))):
                args.append(['python', main_py_path, 
                             os.path.join(dir, str(i)), # samples path
                             str(idx), # unique identifier
                             seg, # mask threshold
                             '1', # excised root
                             '1', # crown root
                             '1', # segmentation
                             '0.0', # marker diameter
                             '0', # stem reconstruction
                             '0', # plots
                             '0', # output format
                             work_dir, # working directory
                             './traits.csv']) # trait file path
    r = pool.map_async(calculate, args)
    r.wait() # Wait on the results
    print 'All files done in '+str(time.time()-startT)+'s !'
    print 'Collecting results'
    oc.combineOutput(work_dir) 
    print 'Results written to ' + os.path.join(work_dir, 'outputAll.csv')
