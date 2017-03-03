#!/usr/bin/env python
import os,sys
from shutil import move
import pyfits
from pyraf import iraf
import tempfile 
from misc import *
from zcombine import * 
from mkobjmask import * 
from optparse import OptionParser

def mkskyflat(inlist, output, inpref='', mskpref='m', bpm='none', iternum=50, nsigma=3, minpix=250, hsigma=3, lsigma=10, conv='block 3 3'):
    
    # load IRAF package
    iraf.proto()
    iraf.nproto()

    # open input list and check if it exists
    inimg_arr = check_input(inlist, inpref)
    if len(inimg_arr) == 1 and inimg_arr[0] == 1:
        return 1

    # check if output image is already exist
    if os.access(output, os.R_OK) : 
        print >> sys.stderr, 'operation would overwrite existing file (%s)' % output
        return 1

    # check optional parameters 
    if bpm.lower() == 'none' or bpm == '':
        bpmex = 0
    else:
        bpmex = 1
        if not os.access(bpm, os.R_OK):
            print >> sys.stderr, 'cannot read bad pixel mask (%s)' % bpm
            return 1

    if iternum <= 0:
        print >> sys.stderr, 'invalid iternum value (%d)' % iternum
        return 1

    if nsigma <= 0:
        print >> sys.stderr, 'invalid nsigma value (%f)' % nsigma
        return 1

    # prefix for temoprary images
    tmp = tempfile.NamedTemporaryFile(suffix='', prefix='', dir='/tmp')
    tmp_prefix = tmp.name
    tmp.close()
    
    # temporary normalize raw frames by average
    tmp_list = tmp_prefix+'.lst'
    fout = open(tmp_list, 'w')
    for i in range(len(inimg_arr)):

        tmpimg = tmp_prefix+os.path.basename(inimg_arr[i])
        fout.write('%s\n' % tmpimg)

        # remove temporary files if it exists
        if os.access(tmpimg, os.R_OK):
            os.remove(tmpimg)

        # calculate average counts in each frame
        if bpmex == 0:
            iraf.unlearn('imstat')
            res = iraf.imstat(inimg_arr[i], format='no', fields='mean', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
            if len(res) == 1:
                ave = float(res[0])
            else:
                print >> sys.stderr, 'failed to calculate avarage from the frame (%s)' % inimg
                remove_temp_all(tmp_prefix)
                fout.close()
                return 1
        else:
            iraf.unlearn('mimstat')
            res = iraf.mimstat(inimg_arr[i], imasks=bpm, format='no', fields='mean', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
            if len(res) == 1:
                ave = float(res[0])
            else:
                print >> sys.stderr, 'failed to calculate masked avarage from the frame (%s)'% inimg
                remove_temp_all(tmp_prefix)
                fout.close()
                return 1

        # normalize each frame by average
        iraf.unlearn('imarith')
        iraf.imarith(inimg_arr[i], '/', ave, tmpimg, verbose='no')

    # close output file handler 
    fout.close()

    # make temporary sky flat frame
    tmp_flat = tmp_prefix+'_sflat.fits'
    ret = zcombine(tmp_list, tmp_flat, combine='median', reject='sigclip')
    if ret != 0:
        remove_temp_all(tmp_prefix)
        return 1

    # temporary flat field 
    iraf.unlearn('imarith')
    tmpff_list = tmp_prefix+'_ff.lst'
    fout = open(tmpff_list, 'w')
    for i in range(len(inimg_arr)):

        tmpimg = tmp_prefix+'_ff_'+os.path.basename(inimg_arr[i])
        fout.write('%s\n' % tmpimg)

        # remove temporary files if it exists
        if os.access(tmpimg, os.R_OK):
            os.remove(tmpimg)

        # divided by temporary flat 
        iraf.imarith(inimg_arr[i], '/', tmp_flat, tmpimg, verbose='no')

        # read image fits header
        im = pyfits.open(inimg_arr[i])
        exptime = float(im[0].header['EXP1TIME'])
        coadds = int(im[0].header['COADDS'])
        ndr = int(im[0].header['NDR'])
        im.close()

        # normalize by coadds, ndr, and exptime values
        iraf.imarith(tmpimg, '/', (coadds*ndr*exptime), tmpimg, verbose='no')


    # close output file handler 
    fout.close()

    # remove temporary flat
    os.remove(tmp_flat)

    # calculate standard deviation of the temporary flat divided images
    try:
        iraf.unlearn('imstat')
        stddev = iraf.imstat('@'+tmpff_list, format='no', fields='stddev', nclip=50, Stdout=1)
    except:
        remove_temp_all(tmp_prefix)
        return 1

    # calculate mean and sigma of the stddev of the image
    y = map(float, stddev)
    mean, sigma = meanclip(np.array(y, dtype=float))
    #print mean, sigma 

    # flag bad frames
    frame_flag = []
    for i in range(len(stddev)):
        im = pyfits.open(inimg_arr[i], mode='update')                
        if math.fabs(y[i] - mean) < nsigma * sigma:
            frame_flag.append(0)
            im[0].header['PATTERN'] = 0
        else:
            frame_flag.append(1)
            im[0].header['PATTERN'] = 1
        im.close()

    # make object mask 
    ret = mkobjmask(tmpff_list, inpref='', outpref=mskpref, bpm=bpm, minpix=minpix, hsigma=hsigma, lsigma=lsigma, conv=conv)

    # rename mask file 
    mskimg_arr = []
    for i in range(len(inimg_arr)):
        base = os.path.basename(inimg_arr[i])
        tmpimg = tmp_prefix+'_ff_'+base
        maskimg = mskpref+base
        tmpmask = mskpref+os.path.basename(tmpimg)
        move(tmpmask, maskimg)
        mskimg_arr.append(maskimg)

    # delete temporary flat fielded files 
    remove_temp_all(tmp_prefix)
    
    # check mkobjmask result 
    if ret != 0:
        return 1

    # nomalize by object masked average 
    fout = open(tmp_list, 'w')
    for i in range(len(inimg_arr)):

        if frame_flag[i] == 0:
        
            tmpimg = tmp_prefix+os.path.basename(inimg_arr[i])
            maskimg = mskimg_arr[i]
            tmpmask = tmp_prefix+maskimg
            fout.write('%s\n' % tmpimg)

            if os.access(tmpimg, os.R_OK):
                os.remove(tmpimg)
                
            if os.access(tmpmask, os.R_OK):
                os.remove(tmpmask)

            if bpmex == 1:
                iraf.unlearn('imarith')
                iraf.imarith(maskimg+'[pl]','+',bpm,tmpmask,verbose='no')
            else :
                iraf.unlearn('imcopy')
                iraf.imcopy(maskimg+'[pl]',tmpmask,verbose='no')

            iraf.unlearn('mimstat')
        
            # calculate average of the raw images with mask
            res = iraf.mimstat(inimg_arr[i], imasks=tmpmask, format='no', fields='mean', nclip=iternum, lsigma=nsigma, usigma=nsigma, Stdout=1)
            if len(res) == 1:
                ave = float(res[0])
            else:
                print >> sys.stderr, 'failed to calculate masked avarage from the frame (%s)' % inimg_arr[i]
                remove_temp_all(tmp_prefix)
                return 1
        
            # normalize each frame by average
            iraf.unlearn('imarith')
            iraf.imarith(inimg_arr[i], '/', ave, tmpimg, verbose='no')
            
            # add mask into header 
            iraf.unlearn('hedit')
            iraf.hedit(tmpimg, 'objmask', maskimg+'[pl]', add='yes', verify='no')

    # close output file handler 
    fout.close()
        
    # combine
    ret = zcombine(tmp_list, output, combine='median', reject='sigclip', masktype='!objmask')
    if ret != 0:
        print >> sys.stderr, 'failed to combine normalized frames'
        remove_temp_all(tmp_prefix)
        return 1   

    # remove temoprary files
    remove_temp_all(tmp_prefix)

    return 0


if __name__=="__main__":

    usage = "usage: %prog input_list output_image [options]"

    parser = OptionParser(usage)
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=null)")
    parser.add_option("--mskpref", dest="mskpref", type="string", default="m",
                      help="prefix for output mask frame or optuput mask list with @ mark (default=m)")
    parser.add_option("--bpm", dest="bpm", type="string", default="none",
                      help="bad pixel mask frame (default=none)")
    parser.add_option("--iternum", dest="iternum", type="int", default=50,
                      help="Number of clipping iterations (default=50)")
    parser.add_option("--nsigma", dest="nsigma", type="float", default=3.0,
                      help="N-sigma reject limit (default=3.0)")
    parser.add_option("--minpix", dest="minpix", type="int", default=250,
                      help="minimum number of pixels in detected objects (default=250)")
    parser.add_option("--hsigma", dest="hsigma", type="float", default=3.0,
                      help="sigma threshold above sky (default=3)")
    parser.add_option("--lsigma", dest="lsigma", type="float", default=10.0,
                      help="sigma threshold below sky (default=10)")
    parser.add_option("--conv", dest="conv", type="string", default="block 3 3",
                      help="convolution kernel for objmask (default=\"block 3 3\")")
    
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()

    mkskyflat(args[0], args[1], inpref=options.inpref, mskpref=options.mskpref, bpm=options.bpm, iternum=options.iternum, nsigma=options.nsigma, minpix=options.minpix, hsigma=options.hsigma, lsigma=options.lsigma, conv=options.conv)
