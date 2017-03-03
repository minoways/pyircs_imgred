#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
import tempfile 
from misc import * 
from optparse import OptionParser
import numpy as np
from numpy import linalg as la

def getLGSpos(inimg):
    
    ### constants ###
    #ttmtx = np.matrix('0.33278 0; 0 -0.33278')
    wttmtx = np.matrix('0.523603 -0.512550; -0.512550 -0.523603')
    #tt0 = np.matrix('0.0; 0.0')
    wtt0 = np.matrix('5.447; 5.004')
    ircs52mtx = np.matrix('-0.028168 0.000058; 0.000059 0.028168')

    # check image
    if os.access(inimg, os.R_OK):
        im = pyfits.open(inimg)
    else:
        print >> sys.stderr, 'cannot read image (%s)' % inimg
        f.close()
        return 1

    # read image data and fits header
    hdr0 = im[0].header
    au1gsx = float(im[0].header['D_AU1GSX'])
    au1gsy = float(im[0].header['D_AU1GSY'])
    wtt1 = float(im[0].header['D_WTTC1'])
    wtt2 = float(im[0].header['D_WTTC2'])
    
    wtt = np.matrix((wtt1, wtt2)).T

    # close input image 
    im.close()

    # calculate WTT offset due to tracking error/LLT error
    dau_wtt_offset = -1.0*wttmtx.dot((wtt-wtt0))
    
    # calculate position offset on the science camera from the AU position
    offset = la.inv(ircs52mtx).dot(dau_wtt_offset)
    
    # calculate LGS position
    xLGS = au1gsx + offset[0]
    yLGS = au1gsy + offset[1]

    print "LGS position (x.y) = ( %.2f , %.2f )" % (xLGS,yLGS)

    return 0

            
if __name__=="__main__":

    usage = "usage: %prog input_image"

    parser = OptionParser(usage)
    
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit()

    getLGSpos(args[0])
   
