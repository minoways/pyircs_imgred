#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
import tempfile 
from misc import *
from optparse import OptionParser

def mkobjmask(inlist, inpref='', outpref='m', bpm='none', minpix=250, hsigma=3, lsigma=10, conv='block 3 3'):
    
    # load IRAF package
    iraf.nproto()

    # open input list and check if it exists
    inimg_arr = check_input(inlist, inpref)
    if isinstance(inimg_arr, int):
        return 1

    # check output images 
    outimg_arr = check_outpref(outpref, inimg_arr)
    if isinstance(outimg_arr, int):
        return 1 

    # check optional parameters 
    if bpm.lower() != 'none':
        if not os.access(bpm, os.R_OK):
            print >> sys.stderr, 'cannot read bad pixel mask (%s)' % bpm
            return 1
        bpmex = 1
    else:
        bpmex = 0
     
    if lsigma <= 0:
        print >> sy.stderr, 'invalid lsigma value (%f)' % lsigma
        return 1

    if hsigma <= 0:
        print >> sys.stderr, 'invalid hsigma value (%f)' % hsigma
        return 1

    # check convolution kernel name
    ret = 1
    param = conv.split()
    if len(param) == 1:
        if os.access(param[0], os.R_OK):
            conv = '@'+conv
            ret = 0
        else:
            print >> sys.stderr, 'convolution kernel file (%s) does not exist' % param[0]
    elif param[0] == 'block' and len(param) == 3:
        if param[1].isdigit() and param[2].isdigit():
            ret = 0
    elif param[0] == 'bilinear' and len(param) == 3:
        if param[1].isdigit() and param[2].isdigit() and int(param[1]) > 0 and int(param[2]) > 0:
            ret = 0
    elif param[0] == 'gauss' and len(param) == 5:
        if param[1].isdigit() and param[2].isdigit() and int(param[1]) > 0 and int(param[2]) > 0:
            if isfloat(param[3]) and isfloat(param[4]) and float(param[3]) > 0 and float(param[4])>0:
                ret = 0

    if ret > 0:
        print >> sys.stderr, 'invalid convolve parameter (%s)' % conv
        return 1

    # prefix for temoprary images
    tmp = tempfile.NamedTemporaryFile(suffix='', prefix='', dir='/tmp')
    tmp_prefix = tmp.name
    tmp.close()
    
    # check if input images exist and fix bad pixel if requested
    for i in range(len(inimg_arr)):
            
        tmp_inimg = tmp_prefix+os.path.basename(inimg_arr[i])
        tmp_outimg = tmp_prefix+'mask'+os.path.basename(outimg_arr[i])
        iraf.unlearn('imcopy')
        iraf.imcopy(inimg_arr[i], tmp_inimg,verbose='no')
        if bpmex == 1:
            iraf.unlearn('fixpix')
            iraf.fixpix(tmp_inimg, bpm)
            iraf.fixpix(tmp_inimg, bpm)
            iraf.fixpix(tmp_inimg, bpm)
        iraf.unlearn('objmasks')
        iraf.objmasks(tmp_inimg, tmp_outimg, omtype='boolean', masks='', convolve=conv, lsigma=lsigma, hsigma=hsigma, hdetect='yes', ldetect='yes', minpix=minpix)
        ret = []
        ret = iraf.objmasks(tmp_inimg, outimg_arr[i], omtype='boolean', masks=tmp_outimg+'[pl]', convolve=conv, lsigma=lsigma, hsigma=hsigma, hdetect='yes', ldetect='yes', minpix=minpix, Stdout=1)
        if len(ret) != 0:
            iraf.unlearn('hedit')
            iraf.hedit(inimg_arr[i], 'objmask', outimg_arr[i]+'[pl]', add='yes', verify='no')
        os.remove(tmp_inimg)
        os.remove(tmp_outimg)
        if len(ret) == 0:
            print >> sys.stderr, 'failed to execute iraf objmasks task for %s' % inimg_arr[i]
            return 1
    
    return 0

if __name__=="__main__":

    usage = "usage: %prog input_list [options]"

    parser = OptionParser(usage)
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=none)")
    parser.add_option("--outpref", dest="outpref", type="string", default="m",
                      help="prefix for output frame or output list with @ mark (default=m)")
    parser.add_option("--bpm", dest="bpm", type="string", default="none",
                      help="bad pixel mask frame (default=none)")
    parser.add_option("--minpix", dest="minpix", type="int", default=250,
                      help="minimum number of pixels in detected objects (default=250)")
    parser.add_option("--hsigma", dest="hsigma", type="float", default=3.0,
                      help="sigma threshold above sky (default=3)")
    parser.add_option("--lsigma", dest="lsigma", type="float", default=10.0,
                      help="sigma threshold below sky (default=10)")
    parser.add_option("--conv", dest="conv", type="string", default="block 3 3",
                      help="convolution kernel for objmask (default=\"block 3 3\")")
    
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit()

    mkobjmask(args[0], inpref=options.inpref, outpref=options.outpref, bpm=options.bpm, minpix=options.minpix, hsigma=options.hsigma, lsigma=options.lsigma, conv=options.conv)
