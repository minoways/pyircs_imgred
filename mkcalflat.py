#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
import tempfile 
from optparse import OptionParser
from misc import *

def mkcalflat(calonlst, calofflst, outimg, mask='none'):
    
    # load IRAF package
    iraf.proto()

    # open cal off list and check if files exist
    on_list = check_input(calonlst, '')
    if len(on_list) == 1 and on_list[0] == 1:
        return 1

    # open cal off list and check if files exist
    off_list = check_input(calofflst, '')
    if len(off_list) == 1 and off_list[0] == 1:
        return 1

    # check if output file exists 
    if os.access(outimg, os.R_OK):
        print >> sys.stderr, 'operation would overwrite existing image (%s)' % outimg
        return 1

    # prefix for temporary images and files
    tmp = tempfile.NamedTemporaryFile(suffix='', prefix='', dir='/tmp')
    tmp_prefix = tmp.name
    tmp.close()

    # check exptime, coadds, ndr, slwcnt for off frames
    tmp_off_list = tmp_prefix + '_off.lst'
    fout = open(tmp_off_list, 'w')
    err = 0
    for i in range(len(off_list)):
        im = pyfits.open(off_list[i])                
        im_exptime = float(im[0].header['EXP1TIME'])
        im_coadds = int(im[0].header['COADDS'])
        im_ndr = int(im[0].header['NDR'])
        im_slwcnt = int(im[0].header['I_SLWCNT'])
        if i == 0:
            off_exptime = im_exptime
            off_coadds = im_coadds
            off_ndr = im_ndr
            off_slwcnt = im_slwcnt
        else:
            if im_exptime != off_exptime:
                print >> sys.stderr, 'exptime for %s is not consistent with the first frame (%s)' % (off_list[i],off_list[0])
                err += 1
            if im_coadds != off_coadds:
                print >> sys.stderr, 'coadds for %s is not consistent with the first frame (%s)' % (off_list[i],off_list[0])
                err += 1
            if im_ndr != off_ndr:
                print >> sys.stderr, 'ndr for %s is not consistent with the first frame (%s)' % (off_list[i],off_list[0])
                err += 1
            if im_slwcnt != off_slwcnt:
                print >> sys.stderr, 'slwcnt for %s is not consistent with the first frame (%s)' % (off_list[i],off_list[0])
                err += 1
        
        im.close()
        fout.write('%s\n' % off_list[i])
    fout.close()

    # check error
    if err != 0:
        remove_temp_all(tmp_prefix)
        return 1

    # check exptime, coadds, ndr, slwcnt for on frames
    tmp_on_list = tmp_prefix + '_on.lst'
    fout = open(tmp_on_list, 'w')
    for i in range(len(on_list)):
        im = pyfits.open(on_list[i])                
        im_exptime = float(im[0].header['EXP1TIME'])
        im_coadds = int(im[0].header['COADDS'])
        im_ndr = int(im[0].header['NDR'])
        im_slwcnt = int(im[0].header['I_SLWCNT'])
        if i == 0:
            on_exptime = im_exptime
            on_coadds = im_coadds
            on_ndr = im_ndr
            on_slwcnt = im_slwcnt
        else:
            if im_exptime != on_exptime:
                print >> sys.stderr, 'exptime for %s is not consistent with the first frame (%s)' % (on_list[i],on_list[0])
                err += 1
            if im_coadds != on_coadds:
                print >> sys.stderr, 'coadds for %s is not consistent with the first frame (%s)' % (on_list[i],on_list[0])
                err += 1
            if im_ndr != on_ndr:
                print >> sys.stderr, 'ndr for %s is not consistent with the first frame (%s)' % (on_list[i],on_list[0])
                err += 1
            if im_slwcnt != on_slwcnt:
                print >> sys.stderr, 'slwcnt for %s is not consistent with the first frame (%s)' % (on_list[i],on_list[0])
                err += 1

        im.close()
        fout.write('%s\n' % on_list[i])
    fout.close()
    
    # check error
    if err != 0:
        remove_temp_all(tmp_prefix)
        return 1

    # check consistency of EXPTIME, COADDS, NDR and SLWCNT between on and off frames
    if on_exptime != off_exptime:
        print >> sys.stderr, 'exptime for on and off frames are not consistent'
        err += 1
    if on_coadds != off_coadds:
        print >> sys.stderr, 'coadds for on and off frames are not consistent'
        err += 1
    if on_ndr != off_ndr:
        print >> sys.stderr, 'ndr for on and off frames are not consistent'
        err += 1
    if on_slwcnt != off_slwcnt:
        print >> sys.stderr, 'slwcnt for on and off frames are not consistent'
        err += 1

    # check error
    if err != 0:
        remove_temp_all(tmp_prefix)
        return 1

    # combine off frames 
    tmpoff = tmp_prefix + '_off.fits'
    iraf.unlearn('imcombine')
    iraf.imcombine.combine = 'median'
    iraf.imcombine.reject = 'sigclip'
    iraf.imcombine('@'+tmp_off_list, tmpoff)

    # combine on frames 
    tmpon = tmp_prefix + '_on.fits'
    tmp.close()
    iraf.unlearn('imcombine')
    iraf.imcombine.combine = 'median'
    iraf.imcombine.reject = 'sigclip'
    iraf.imcombine('@'+tmp_on_list, tmpon)
    
    # subtract off from on
    iraf.unlearn('imarith')
    iraf.imarith(tmpon,"-",tmpoff,outimg,verbose='no')

    # delete temporary files 
    remove_temp_all(tmp_prefix)

    # calculate average 
    if mask == 'none':
        iraf.unlearn('imstat')
        res = iraf.imstat(outimg, format='no', fields='mean', nclip=50, Stdout=1)
        if len(res) == 1:
            ave = float(res[0])
        else:
            print >> sys.stderr, 'failed to calculate avarage of the combined frame'
            os.remove(outimg)
            return 1
    else:
        if os.access(mask, os.R_OK):
            iraf.unlearn('mimstat')
            res = iraf.mimstat(outimg, imasks=mask, format='no', fields='mean', nclip=50, Stdout=1)
            if len(res) == 1:
                ave = float(res[0])
            else:
                print >> sys.stderr, 'failed to calculate masked avarage of the combined frame' 
                os.remove(outimg)
                return 1
        else :
            print >> sys.stderr, 'mask frame (%s) does not exist' % mask
            os.remove(outimg)
            return 1 

    # divide by average count
    iraf.unlearn('imarith')
    iraf.imarith(outimg,"/",ave,outimg,verbose='no')

    return 0

            
if __name__=="__main__":

    usage = "usage: %prog cal_on_list cal_off_list output_img [options]"

    parser = OptionParser(usage)
    parser.add_option("--bpm", dest="bpm", type="string", default="none",
                      help="bad pixel mask frame (default=none)")
    
    (options, args) = parser.parse_args()
    if len(args) != 3:
        parser.print_help()
        sys.exit()

    mkcalflat(args[0],args[1],args[2],mask=options.bpm)
   
