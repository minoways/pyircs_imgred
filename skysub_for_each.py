#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
import tempfile 
from zcombine import * 
from misc import *
from operator import itemgetter 
from optparse import OptionParser

def skysub_for_each(inimg, mskinimg, skylist, outpref='s', mskpref='none', iternum=50, nsigma=3, second=False, step=0.05):
    
    # load IRAF package
    iraf.proto()

    # check input image
    if not os.access(inimg, os.R_OK):
        print >> sys.stderr, 'cannot read input image (%s)' % inimg
        return 1

    # check mask image for input 
    if mskinimg.lower() == 'none' or mskinimg == '':
        mask_inimg = 0
    else:
        mask_inimg = 1
        if not os.access(mskinimg, os.R_OK):
            print >> sys.stderr, 'cannot read input mask image (%s)' % mskinimg
            return 1

    # check output
    outimg = outpref + inimg
    if os.access(outimg, os.R_OK):
        print >> sys.stderr, 'operation would overwrite existing image (%s)' % outimg
        return 1

    # check input sky image list 
    skyimg_arr = check_input(skylist, '')
    if isinstance(skyimg_arr, int):
        return 1

    # check mask images for sky frames 
    if mskpref.lower() == 'none' or mskpref == '':
        mask = 0
    else:
        mask = 1
        mskimg_arr = check_inpref(mskpref, skyimg_arr)
        if isinstance(mskimg_arr, int):
            return 1

    # check step size
    if step >= 0.5:
        print >> sys.stderr, 'step size should be smaller than 0.5'
        return 1

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
    for i in range(len(skyimg_arr)):

        if mask == 1:
            if mskpref[0] != ':':
                if not second:
                    iraf.hedit(inimg_arr[i], 'objmask', mskimg_arr[i]+'[pl]', add='yes', verify='no')
                else:
                    iraf.hedit(inimg_arr[i], 'objmask', mskimg_arr[i], add='yes', verify='no')



    # sky subtraction 
    tmp = tempfile.NamedTemporaryFile(suffix='', prefix='', dir='/tmp')
    tmp_prefix = tmp.name
    tmp.close()

    tmp_sky = '%s_sky_%03d.fits' % (tmp_prefix, int(i))
    # combine subset 
    if mask == 1:
        ret = zcombine(skylist, tmp_sky, combine='median', reject='sigclip', masktype='!OBJMASK')
    else:
        ret = zcombine(skylist, tmp_sky, combine='median', reject='sigclip', masktype='none')
    if ret != 0:
        remove_temp_all(tmp_prefix)
        return 1
    
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
    sky_arr = np.reshape(np.array(sky[0].data - skymed), sky[0].data.size)
    sky.close()
    os.remove(tmp_sky)

    # subtract sky
    if mask_inimg == 1:
        iraf.unlearn('mimstat')
        if not second:
            ret = iraf.mimstat(inimg, imasks=mskinimg+'[pl]', format='no', fields='midpt', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
        else:
            ret = iraf.mimstat(inimg, imasks=mskinimg, format='no', fields='midpt', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
    else: 
        iraf.unlearn('imstat')
        ret = iraf.imstat(inimg, format='no', fields='midpt', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
    
    if len(ret) == 1:
        objmed = float(ret[0])
    else : 
        remove_temp_all(tmp_prefix)
        print >> sys.stderr, 'failed to calculate median of the object background'
        return 1

    # record input fits data into memory
    im = pyfits.open(inimg)
    im_hdr = im[0].header
    im_data = im[0].data
    im_arr = np.reshape(np.array(im[0].data - objmed), im[0].data.size)
    im.close()
            
    # check if minimum sigma is obtained at scale_factor = 1.0
    m1,s1 = meanclip(im_arr - sky_arr)
    m2,s2 = meanclip(im_arr - (1.0+step)*sky_arr)
    m3,s3 = meanclip(im_arr - (1.0-step)*sky_arr)

    if s1 < s2 and s1 < s3:
        print 'scale factor = 1.0'
        # subtract sky in memory 
        out_data = (im_data - objmed) - (sky_data - skymed)

    else:
        # find minimum sigma with coarse sampling 
        n = 1
        min_i = 0
        while min_i - 1 < 0 or min_i + 1 >= len(s_list):
            s = -1.0 * n 
            ds = 0.5
            sigma_list=[]
            s_list=[]
            while s <= 1.0*n:
                sub = im_arr - (s * sky_arr)
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
            sub = im_arr - (s * sky_arr)
            mean,sigma=meanclip(sub)
            sigma2_list.append(sigma)
            s2_list.append(s)
            s += ds
        min_i, min_sigma= min(enumerate(sigma2_list),key=itemgetter(1))
    
        print 'scale factor = %f' % s2_list[min_i]
            
        # subtract sky in memory 
        out_data = (im_data - objmed) - (sky_data - skymed)*s2_list[min_i]
            
    # save output fits file 
    hdu = pyfits.PrimaryHDU(out_data)
    imgout = pyfits.HDUList([hdu])
    imgout[0].header = im_hdr
    imgout.writeto(outimg)
    imgout.close()

    remove_temp_all(tmp_prefix)
    return 0

if __name__=="__main__":

    usage = "usage: %prog input_image input_mask sky_list [options]"

    parser = OptionParser(usage)
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
    parser.add_option("--step", dest="step", type="float", default=0.05,
                      help="step for finding the best scale value (default=0.05)")

    (options, args) = parser.parse_args()
    if len(args) != 3:
        parser.print_help()
        sys.exit()
                
    skysub_for_each(args[0], args[1], args[2], outpref=options.outpref, mskpref=options.mskpref, iternum=options.iternum, nsigma=options.nsigma, second=options.second, step=options.step)
  
