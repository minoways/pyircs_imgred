#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
import tempfile 
from zcombine import * 
from misc import *
from optparse import OptionParser

def skysub(inlist, inpref='', outpref='s', mskpref='none', ndthpos=9, iternum=50, nsigma=3, second=False):
    
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

    # check dithering width
    if ndthpos != 5 and ndthpos != 9:
        print >> sys.stderr, 'number of dithering position should be 5 or 9' 

    if iternum <= 0:
        print >> sys.stderr, 'invalid iternum value (%d)' % iternum
        return 1

    if nsigma <= 0:
        print >> sys.stderr, 'invalid nsigma value (%f)' % nsigma
        return 1

    # check input, output, and mask files and record mask frame into header
    pre_pos = 100
    in_array = []
    out_array = []
    mask_array = []
    iraf.unlearn('hedit')
    for i in range(len(inimg_arr)):

        if mask == 1:
            if mskpref[0] != ':':
                if not second:
                    iraf.hedit(inimg_arr[i], 'objmask', mskimg_arr[i]+'[pl]', add='yes', verify='no')
                else:
                    iraf.hedit(inimg_arr[i], 'objmask', mskimg_arr[i], add='yes', verify='no')

        # get dithering position from header 
        im = pyfits.open(inimg_arr[i])
        pos = int(im[0].header['I_DTHPOS'])
        im.close()

        # record each dithering set into subset 
        if pos < pre_pos:
            if i != 0:
                in_array.append(in_subset)
                out_array.append(out_subset)
                if mask == 1:
                    mask_array.append(mask_subset)
            in_subset = []
            out_subset = []
            in_subset.append(inimg_arr[i])
            out_subset.append(outimg_arr[i])
            if mask == 1:
                mask_subset = []
                mask_subset.append(mskimg_arr[i])
            pre_pos = pos
        else:
            in_subset.append(inimg_arr[i])
            out_subset.append(outimg_arr[i])
            if mask == 1:
                mask_subset.append(mskimg_arr[i])
            pre_pos = pos
    
    # record rest of frames
    in_array.append(in_subset)
    out_array.append(out_subset)
    if mask == 1:
        mask_array.append(mask_subset)

    # number of subsets 
    n_subset = len(in_array)

    # sky subtraction 
    tmp = tempfile.NamedTemporaryFile(suffix='', prefix='', dir='/tmp')
    tmp_prefix = tmp.name
    tmp.close()

    for i in range(n_subset):
        tmp_sky = '%s_sky_%03d.fits' % (tmp_prefix, int(i))
        tmp_subset = '%s_subset_%03d.lst' % (tmp_prefix, int(i))
        # make subset list 
        f = open(tmp_subset, 'w')
        for j in range(len(in_array[i])):
            f.write(in_array[i][j]+'\n')
        f.close()
        # combine subset 
        if mask == 1:
            ret = zcombine(tmp_subset, tmp_sky, combine='median', reject='sigclip', masktype='!OBJMASK')
        else:
            ret = zcombine(tmp_subset, tmp_sky, combine='median', reject='sigclip', masktype='none')
        if ret != 0:
            remove_temp_all(tmp_prefix)
            return 1
        os.remove(tmp_subset)

        # calculate median of combined sky background 
        iraf.unlearn('imstat')
        ret = iraf.imstat(tmp_sky, format='no', fields='midpt', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
        if len(ret) == 1:
            skymed = float(ret[0])
        else : 
            remove_temp_all(tmp_prefix)
            print >> sys.stderr, 'failed to calculate median of the combined sky background'
            return 1

        # save sky data into memory
        sky = pyfits.open(tmp_sky)
        sky_data = sky[0].data
        sky.close()
        os.remove(tmp_sky)

        # subtract sky
        for j in range(len(in_array[i])):
            if mask == 1:
                iraf.unlearn('mimstat')
                if not second:
                    ret = iraf.mimstat(in_array[i][j], imasks=mask_array[i][j]+'[pl]', format='no', fields='midpt', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
                else:
                    ret = iraf.mimstat(in_array[i][j], imasks=mask_array[i][j], format='no', fields='midpt', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
            else: 
                iraf.unlearn('imstat')
                ret = iraf.imstat(in_array[i][j], format='no', fields='midpt', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
            if len(ret) == 1:
                objmed = float(ret[0])
            else : 
                remove_temp_all(tmp_prefix)
                print >> sys.stderr, 'failed to calculate median of the object background'
                return 1
            # record input fits data into memory
            im = pyfits.open(in_array[i][j])
            im_hdr = im[0].header
            im_data = im[0].data
            im.close()
            
            # subtract sky in memory 
            out_data = im_data - (sky_data * objmed / skymed)
            
            # save output fits file 
            hdu = pyfits.PrimaryHDU(out_data)
            imgout = pyfits.HDUList([hdu])
            imgout[0].header = im_hdr
            imgout.writeto(out_array[i][j])
            imgout.close()

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
    parser.add_option("--ndthpos", dest="ndthpos", type="int", default=9,
                      help="number of dithering position (default=9)")
    parser.add_option("--iternum", dest="iternum", type="int", default=50,
                      help="Number of clipping iterations (default=50)")
    parser.add_option("--nsigma", dest="nsigma", type="float", default=3.0,
                      help="N-sigma reject limit (default=3.0)")
    parser.add_option("--second", dest="second", action="store_true", default=False,
                      help="Is this second cycle subtraction? (default=False)")
    
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit()
                
    skysub(args[0], inpref=options.inpref, outpref=options.outpref, mskpref=options.mskpref, ndthpos=options.ndthpos, iternum=options.iternum, nsigma=options.nsigma, second=options.second)
  
