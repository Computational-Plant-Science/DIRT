'''
Masking.py

The MAsking module for DIRT. This class contains all functions to compute the binary masked as described in the paper.

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
 external library imports
'''
import mahotas as m
import numpy as np
import scipy.ndimage

class Masking(object):
    '''
    classdocs
    '''


    def __init__(self,scale=1.0):
        '''
        Constructor
        '''
        self.__scale=scale
        

    def threshold_adaptive(self,image, block_size, method='gaussian', offset=0,
                       mode='reflect', param=None):
        thresh_image = np.zeros(image.shape, 'double')
        if method == 'generic':
            scipy.ndimage.generic_filter(image, param, block_size,
                output=thresh_image, mode=mode)
        elif method == 'gaussian':
            if param is None:
                # automatically determine sigma which covers > 99% of distribution
                sigma = (block_size - 1) / 6.0
            else:
                sigma = param
            scipy.ndimage.gaussian_filter(image, sigma, output=thresh_image,
                mode=mode)
        elif method == 'mean':
            mask = 1. / block_size * np.ones((block_size,))
            # separation of filters to speedup convolution
            scipy.ndimage.convolve1d(image, mask, axis=0, output=thresh_image,
                mode=mode)
            scipy.ndimage.convolve1d(thresh_image, mask, axis=1,
                output=thresh_image, mode=mode)
        elif method == 'median':
            scipy.ndimage.median_filter(image, block_size, output=thresh_image,
                mode=mode)
    
        return image > (thresh_image - offset)
    
    def calculateMask(self,img):
        print 'Masking input'
        if len(np.unique(img))<=2:
            print 'Binary input detected, no thresholding performed'
            idx1=np.where(img==np.unique(img)[0])
            idx2=np.where(img==np.unique(img)[1])
            img[idx1]=False
            img[idx2]=True
        else:
            print 'Grey input detected'
            T=m.otsu(img,ignore_zeros=False)
            T=T*self.__scale
            img = self.threshold_adaptive(img, 80, 'gaussian',offset=-20,param=T)
            img = m.morph.open(img)

        img = m.morph.close(img)
        ''' just a quick fix of the dilation function that caused the binary image to consist of 0 and 2. Now It should be a real binary image '''
        idx1=np.where(img==np.unique(img)[0])
        idx2=np.where(img==np.unique(img)[1])
        img[idx1]=0
        img[idx2]=255
        
        w,h=np.shape(img)
        img[0,:]=0
        img[:,0]=0
        img[w-1,:]=0
        img[:,h-1]=0
        return img
     
    
