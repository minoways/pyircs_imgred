#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
import tempfile 
from zcombine import * 
from misc import *
from optparse import OptionParser

def skysub2(inlist, skyfits, inpref='', outpref='s', mskpref='none', iternum=50, nsigma=3, second=False):
    
    # load IRAF package
    iraf.proto()

    # check input image list 
    inimg_arr = check_input(inlist, inpref)
    if isinstance(inimg_arr, int):
        return 1

    # check skyfits file
    if not os.access(skyfits,os.R_OK):
        print >> sys.stderr, 'cannot open sky fits frame (%s)' % skyfits
        return 1

    # check output 
    outimg_arr = check_outpref(outpref, inimg_arr)
    if isinstance(outimg_arr, int):
        return 1

    # check mask images 
    if mskpref.lower() == 'none' or mskpref == '':
        mask = 0
    else:
        mask = 1
        mskimg_arr = check_inpref(mskpref, inimg_arr)
        if isinstance(mskimg_arr, int):
            return 1

    if iternum <= 0:
        print >> sys.stderr, 'invalid iternum value (%d)' % iternum
        return 1

    if nsigma <= 0:
        print >> sys.stderr, 'invalid nsigma value (%f)' % nsigma
        return 1

        
    # calculate median of combined sky background 
    iraf.unlearn('imstat')
    ret = iraf.imstat(skyfits, format='no', fields='midpt', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
    if len(ret) == 1:
        skymed = float(ret[0])
    else : 
        print >> sys.stderr, 'failed to calculate median of the combined sky background'
        return 1

    # save sky data into memory
    sky = pyfits.open(skyfits)
    sky_data = sky[0].data
    sky.close()

    # subtract sky
    for i in range(len(inimg_arr)):
        if mask == 1:
            iraf.unlearn('mimstat')
            if not second:
                ret = iraf.mimstat(inimg_arr[i], imasks=mskimg_arr[i]+'[pl]', format='no', fields='midpt', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
            else:
                ret = iraf.mimstat(inimg_arr[i], imasks=mskimg_arr[i], format='no', fields='midpt', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
        else: 
            iraf.unlearn('imstat')
            ret = iraf.imstat(inimg_arr[i], format='no', fields='midpt', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)

        if len(ret) == 1:
                objmed = float(ret[0])
        else : 
            print >> sys.stderr, 'failed to calculate median of the object background'
            return 1

        # record input fits data into memory
        im = pyfits.open(inimg_arr[i])
        im_hdr = im[0].header
        im_data = im[0].data
        im.close()
            
        # subtract sky in memory 
        out_data = im_data - (sky_data * objmed / skymed)
            
        # save output fits file 
        hdu = pyfits.PrimaryHDU(out_data)
        imgout = pyfits.HDUList([hdu])
        imgout[0].header = im_hdr
        imgout.writeto(outimg_arr[i])
        imgout.close()

    return 0

if __name__=="__main__":

    usage = "usage: %prog input_list skyfits [options]"

    parser = OptionParser(usage)
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=none)")
    parser.add_option("--outpref", dest="outpref", type="string", default="s",
                      help="prefix for output frame or output list with @ mark (default=s)")
    parser.add_option("--mskpref", dest="mskpref", type="string", default="none",
                      help="prefix for mask image, image list with \'@\', or header keyword with \':\'")
    parser.add_option("--iternum", dest="iternum", type="int", default=50,
                      help="Number of clipping iterations (default=50)")
    parser.add_option("--nsigma", dest="nsigma", type="float", default=3.0,
                      help="N-sigma reject limit (default=3.0)")
    parser.add_option("--second", dest="second", action="store_true", default=False,
                      help="Is this second cycle subtraction? (default=False)")
    
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()
                
    skysub2(args[0], args[1], inpref=options.inpref, outpref=options.outpref, mskpref=options.mskpref, iternum=options.iternum, nsigma=options.nsigma, second=options.second)
  
