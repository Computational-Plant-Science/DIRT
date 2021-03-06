'''
Created on Feb 12, 2013

@author: koalaspirit
'''
from pytesser import *
import PIL as pil
import numpy as np
import subprocess
import os
import scipy.misc


def getTextFromImage(tagImg,scratchPath,scratchText='temp'):
    
    scipy.misc.imsave(scratchPath+'temp.bmp',tagImg)
    set_scratch_text_name_root(scratchPath, 'temp.bmp')
    text = image_file_to_string(scratchPath+'temp.bmp', cleanup = cleanup_scratch_flag, graceful_errors=True)
    text = text.translate(None, ",!.;:'{}[]-=()*&^%$#@!~`<>?/|\_+")
    text = ''.join(c for c in text if (c.isalnum() or ' ' or ','))
    text = ' '.join(text.split()) 
    print 'Experiment code: '+text
    return text

def getCodeFromImage(tagImg,scratchPath):
    
    scipy.misc.imsave(scratchPath+'temp.bmp',tagImg)
    args = ['/usr/local/bin/zbarimg','-q',scratchPath+'temp.bmp']
    try:
        code = subprocess.Popen(args, stdout=subprocess.PIPE).communicate()[0]
    except Exception as ex:
        print 'Exception while running zbarimg', ex
    print 'BarCode detected: '+str(code)
    return code
