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

def mkskyflat_with_mask(inlist, output, inpref='', mskpref='m', bpm='none', iternum=50, nsigma=3):
    
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

    # check mask images 
    if mskpref.lower() == 'none' or mskpref == '':
        print >> sys.stderr, 'mask frames are required for this task'
        return 1
    else:
        mskimg_arr = check_inpref(mskpref, inimg_arr)
        if isinstance(mskimg_arr, int):
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
    
    # nomalize by object masked average 
    tmp_list = tmp_prefix+'.lst'
    fout = open(tmp_list, 'w')
    for i in range(len(inimg_arr)):
        
        im = pyfits.open(inimg_arr[i], mode='update')
        try:
            pattern = im[0].header['PATTERN']
        except KeyError:
            pattern = 0
        im.close()

        if pattern == 0:
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
                #iraf.imarith(maskimg+'[pl]','+',bpm,tmpmask,verbose='no')
                iraf.imarith(maskimg,'+',bpm,tmpmask,verbose='no')
            else :
                iraf.unlearn('imcopy')
                #iraf.imcopy(maskimg+'[pl]',tmpmask,verbose='no')
                iraf.imcopy(maskimg,tmpmask,verbose='no')

            iraf.unlearn('mimstat')
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
            if mskpref[0] != ':':
                iraf.hedit(tmpimg, 'objmask', maskimg, add='yes', verify='no')
    
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
                      help="prefix for input mask frame or input mask list with @ mark (default=m)")
    parser.add_option("--bpm", dest="bpm", type="string", default="none",
                      help="bad pixel mask frame (default=none)")
    parser.add_option("--iternum", dest="iternum", type="int", default=50,
                      help="Number of clipping iterations (default=50)")
    parser.add_option("--nsigma", dest="nsigma", type="float", default=3.0,
                      help="N-sigma reject limit (default=3.0)")

    
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()

    mkskyflat_with_mask(args[0], args[1], inpref=options.inpref, mskpref=options.mskpref, bpm=options.bpm, iternum=options.iternum, nsigma=options.nsigma)
