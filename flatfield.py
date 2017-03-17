#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
import tempfile 
from misc import * 
from optparse import OptionParser

def flatfield(inlist, flatimg, inpref='', outpref='f', dark='none', bpm='none'):
    
    # open input list and check if it exists
    inimg_arr = check_input(inlist, inpref)
    if isinstance(inimg_arr, int):
        return 1

    # check flat image
    if os.access(flatimg, os.R_OK):
        im_flat = pyfits.open(flatimg)
        flat_data = im_flat[0].data
        flat_shape = flat_data.shape
        im_flat.close()      
    else:
        print >> sys.stderr, 'cannot read flat image (%s)' % flatimg
        f.close()
        return 1

    # check dark frame 
    if dark.lower() != 'none':
        subdark = 1
        if os.access(dark, os.R_OK):
            im_dark = pyfits.open(dark)
            dark_data = im_dark[0].data
            dark_shape = dark_data.shape
            im_dark.close()
        else : 
            print >> sys.stderr, 'cannot read dark image (%s)' % dark
            return 1
    else:
        subdark = 0

    # check bad pixel mask
    if bpm.lower() != 'none':
        bpmex = 1
        if os.access(bpm, os.R_OK) == False:
            print >> sys.stderr, 'cannot read bad pixel mask (%s)' % bpm
            f.close()
            return 1
    else:
        bpmex = 0

    # check output prefix
    outimg_arr = check_outpref(outpref, inimg_arr)
    if isinstance(outimg_arr, int):
        return 1

    # check if file exist and record exptime, coadds, ndr, slwcnt
    for i in range(len(inimg_arr)):
        
        # read image data and fits header
        im = pyfits.open(inimg_arr[i])                
        hdr0 = im[0].header
        im_shape = im[0].data.shape
        exptime = float(im[0].header['EXP1TIME'])
        coadds = int(im[0].header['COADDS'])
        ndr = int(im[0].header['NDR'])
        
        # normalize by coadds and ndr values
        out_data = im[0].data / (coadds*ndr)

        # close input image 
        im.close()

        # subtract dark if requested 
        if subdark == 1:
            if im_shape[0] == dark_shape[0] and im_shape[1] == dark_shape[1]:
                out_data = out_data - dark_data
            else: 
                print >> sys.stderr, 'input image (%s) and dark (%s) have different image shape' % (inimg, dark)
                return 1
                
        # divide by flat frame
        if im_shape[0] == flat_shape[0] and im_shape[1] == flat_shape[1]:
            out_data = out_data / flat_data
        else: 
            print >> sys.stderr, 'input image (%s) and flat (%s) have different image shape' % (inimg, flatimg)
            return 1

        # normalized by EXP1TIME
        out_data = out_data / exptime

        # write output data into fits 
        hdu = pyfits.PrimaryHDU(out_data)
        imgout = pyfits.HDUList([hdu])
        imgout[0].header = hdr0
        imgout[0].header['EXPTIME'] = 1.0
        imgout[0].header['EXPTORG'] = exptime * coadds
        imgout.writeto(outimg_arr[i])
        imgout.close()

        if bpmex == 1: 
            iraf.unlearn('fixpix')
            iraf.fixpix(outimg_arr[i], bpm, verbose='no')
            iraf.fixpix(outimg_arr[i], bpm, verbose='no')
            iraf.fixpix(outimg_arr[i], bpm, verbose='no')
            
    return 0

            
if __name__=="__main__":
 
    usage = "usage: %prog input_list flat_frame [options]"

    parser = OptionParser(usage)
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=none)")
    parser.add_option("--outpref", dest="outpref", type="string", default="f",
                      help="prefix for output frame or output list with @ mark (default=f)")
    parser.add_option("--dark", dest="dark", type="string", default="none",
                      help="dark frame normalized by coadds and nrd (default=none)")
    parser.add_option("--bpm", dest="bpm", type="string", default="none",
                      help="bad pixel mask frame (default=none)")
    
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()

    flatfield(args[0],args[1], inpref=options.inpref, outpref=options.outpref, dark=options.dark, bpm=options.bpm)
   
