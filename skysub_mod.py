#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
import tempfile 
from zcombine import * 
from misc import *
from operator import itemgetter 
from optparse import OptionParser

def skysub_mod(inlist, inpref='', outpref='s', mskpref='none', ndthpos=9, iternum=50, nsigma=3, second=False, step=0.05):
    
    # load IRAF package
    iraf.proto()

    # check step size
    if step >= 0.5:
        print >> sys.stderr, 'step size should be smaller than 0.5'
        return 1

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
        sky_arr = np.reshape(np.array(sky[0].data - skymed), sky[0].data.size)
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
    parser.add_option("--step", dest="step", type="float", default=0.05,
                      help="step for finding the best scale value (default=0.05)")

    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit()
                
    skysub_mod(args[0], inpref=options.inpref, outpref=options.outpref, mskpref=options.mskpref, ndthpos=options.ndthpos, iternum=options.iternum, nsigma=options.nsigma, second=options.second, step=options.step)
  
