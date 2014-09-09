'''
runOnFolder.py

The module allows to run the DIRT software on a folder with many images in it.

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
    files=os.listdir(dir)
    pool = multiprocessing.Pool(processes=6)
    args=[]
    for idx,i in enumerate(files):
        if i !='.DS_Store':
            if os.path.isfile(dir+i):
                if os.path.isdir(dir+str(idx))==False:
                    args.append(['python', os.getcwd()+'/main.py', dir+'/'+str(i),str(idx),'20.0', '0', '1', '1', '0.0', '0',dir+'/'])
    r = pool.map_async(calculate, args)
    r.wait() # Wait on the results
    print 'All files done in '+str(time.time()-startT)+'s !'
    print 'Collecting results'
    oc.combineOutput(dir)
    print 'Results written to '+dir+'outputAll.csv'
