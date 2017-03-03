#!/usr/bin/env python 

import os
import pyfits
import math
from pyraf import iraf
import tempfile
import numpy as np
from operator import itemgetter 
from misc import *
from optparse import OptionParser

def remove_pattern(inlist, template, inpref='', step=0.05):

    # open input list and check if it exists
    if os.path.splitext(inlist)[-1] == ".fits":
        inimg_arr=[inlist]
    else:
        inimg_arr = check_input(inlist, inpref)
        if isinstance(inimg_arr, int):
            return 1

    # check template image
    if os.path.isfile(template) != True:
        print >> sys.stderr, 'template image (%s) does not exist' % template
        return 1

    # check step size
    if step >= 0.5:
        print >> sys.stderr, 'step size should be smaller than 0.5'
        return 1
    
    # open template image
    try:
        tmp_im = pyfits.open(template)
    except:
        print >> sys.stderr, 'cannot read the template image (%s)' % template
        return 1
    tmp_arr = np.reshape(np.array(tmp_im[0].data), tmp_im[0].data.size)
    tmp_im.close()


    # prefix for temoprary images
    tmp = tempfile.NamedTemporaryFile(suffix='', prefix='', dir='/tmp')
    tmp_prefix = tmp.name
    tmp.close()

    # find scaling value and subtract the scaled template
    for i in range(len(inimg_arr)):
        # open image
        try:
            im = pyfits.open(inimg_arr[i])
        except:
            print >> sys.stderr, 'cannot read the image (%s)' % inimg_arr[i]
            return 1
        im_arr = np.reshape(np.array(im[0].data), im[0].data.size)
        im.close()

        # find minimum sigma with coarse sampling 
        n = 1
        min_i = 0
        while min_i - 1 < 0 or min_i + 1 >= len(s_list):
            s = -1.0 * n 
            ds = 0.5
            sigma_list=[]
            s_list=[]
            while s <= 1.0*n:
                sub = im_arr - (s * tmp_arr)
                mean,sigma=meanclip(sub)
                sigma_list.append(sigma)
                s_list.append(s)
                s += ds
            min_i, min_sigma= min(enumerate(sigma_list),key=itemgetter(1))
            n += 1

        # find minimum sigma with fine sampling
        s = s_list[min_i-1]
        ds = step
        s2_list=[]
        sigma2_list=[]
        while s <= s_list[min_i+1]:
            sub = im_arr - (s * tmp_arr)
            mean,sigma=meanclip(sub)
            sigma2_list.append(sigma)
            s2_list.append(s)
            s += ds
        min_i, min_sigma= min(enumerate(sigma2_list),key=itemgetter(1))
    
        #print scale
        print 'scale factor = %f' % s2_list[min_i]

        # subtract the scaled template
        tmp_scaled = '%s_scaled.fits' % tmp_prefix
        iraf.imarith(template,"*",s2_list[min_i], tmp_scaled, verbose='no')
        iraf.imcopy(inimg_arr[i], inimg_arr[i].replace(".fits","")+"_org.fits", verbose='no')
        iraf.imarith(inimg_arr[i],"-", tmp_scaled, inimg_arr[i], verbose='no')
        
        # remove temporary file
        remove_temp_all(tmp_prefix)
    
    return 0
    
if __name__=="__main__":

    usage = "usage: %prog input_list template [options]"

    parser = OptionParser(usage)
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=none)")
    parser.add_option("--step", dest="step", type="float", default=0.05,
                      help="step for finding the best scale value (default=0.05)")
     
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()

    remove_pattern(args[0], args[1], inpref=options.inpref, step=options.step)
