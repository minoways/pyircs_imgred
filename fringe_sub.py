#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
import tempfile 
from zcombine import * 
from misc import *
from operator import itemgetter 
from optparse import OptionParser

def fringe_sub(inlist, inpref='', outpref='s', mskpref='none', iternum=50, nsigma=3, bpm='none', second=False, fr_x0=0.0, fr_y0=0.0, fr_step=5):
    
    # load IRAF package
    iraf.proto()

    # check input image list 
    inimg_arr = check_input(inlist, inpref)
    if isinstance(inimg_arr, int):
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

    # check bad pixel mask
    if bpm.lower() != 'none':
        bpmex = 1
        if os.access(bpm, os.R_OK) == False:
            print >> sys.stderr, 'cannot read bad pixel mask (%s)' % bpm
            return 1
    else:
        bpmex = 0

    # prefix for temoprary images
    tmp = tempfile.NamedTemporaryFile(suffix='', prefix='', dir='/tmp')
    tmp_prefix = tmp.name
    tmp.close()

    # subtract sky
    for i in range(len(inimg_arr)):
        if mask == 1:
            iraf.unlearn('mimstat')
            if not second:
                ret = iraf.mimstat(inimg_arr[i], imasks=mskimg_arr[i]+'[pl]', format='no', fields='midpt', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
            else:
                if bpmex == 1:
                    tmp_mask = tmp_prefix + os.path.basename(mask_img_arr[i])
                    iraf.imarith(maskimg_arr[i], '+', bpm, tmp_mask, verbose='no') 
                    convert_maskfits_int(tmp_mask, tmp_mask)
                    ret = iraf.mimstat(inimg_arr[i], imasks=tmp_mask, format='no', fields='midpt', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
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
        out_data = im_data - objmed
        im.close()
            
        # fringe subtraction for NB imaging 
        if mask == 1:
            if bpmex == 1:
                tmp_mask = tmp_prefix + os.path.basename(mask_img_arr[i])
                mask_im = pyfits.open(tmp_mask)
            else:
                mask_im = pyfits.open(mskimg_arr[i])
            if not second:
                mask_data = mask_im[1].data
            else:
                mask_data = mask_im[0].data
            mask_im.close()
                
        x = np.ones(out_data.shape[0]*out_data.shape[1]).reshape(out_data.shape[0],out_data.shape[1])* np.arange(out_data.shape[1]) + 1
        y = (np.ones(out_data.shape[1]*out_data.shape[0]).reshape(out_data.shape[1],out_data.shape[0]) * np.arange(out_data.shape[0]) + 1).T
        r = np.sqrt((x-fr_x0)*(x-fr_x0) + (y-fr_y0)*(y-fr_y0))
        rmax = np.max(r)
        rmin = np.min(r)

        r1 = rmin
        while r1 <= rmax:
            r2 = r1 + fr_step
            idx = np.logical_and(r>=r1, r<r2)
            out_subset = out_data[idx]

            if mask == 1:
                mask_subset = mask_data[idx]
                idx_sky = np.where(mask_subset == 0)
                med, sigma = meanclip(out_subset[idx_sky],median=1)
            else:
                med, sigma = meanclip(out_data[idx],median=1)

            out_data[idx] = out_data[idx] - med
            r1 += fr_step

        # save output fits file 
        hdu = pyfits.PrimaryHDU(out_data)
        imgout = pyfits.HDUList([hdu])
        imgout[0].header = im_hdr
        imgout.writeto(outimg_arr[i])
        imgout.close()
        
        if bpmex == 1:
            remove_temp_all(tmp_prefix)

    return 0

if __name__=="__main__":

    usage = "usage: %prog input_list [options]"

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
    parser.add_option("--bpm", dest="bpm", type="string", default="none",
                      help="bad pixel mask frame (default=none)")
    parser.add_option("--second", dest="second", action="store_true", default=False,
                      help="Is this second cycle subtraction? (default=False)")
    parser.add_option("--fr_x0", dest="fr_x0", type="float", default=0.0,
                      help="Fringe center X coordinate (default=0.0)")
    parser.add_option("--fr_y0", dest="fr_y0", type="float", default=0.0,
                      help="Fringe center Y coordinate (default=0.0)")
    parser.add_option("--fr_step", dest="fr_step", type="int", default=5,
                      help="Step size for fringe fitting (default=5)")

    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit()
                
    fringe_sub(args[0], inpref=options.inpref, outpref=options.outpref, mskpref=options.mskpref, iternum=options.iternum, nsigma=options.nsigma, bpm=options.bpm, second=options.second, fr_x0=options.fr_x0, fr_y0=options.fr_y0, fr_step=options.fr_step)
  
