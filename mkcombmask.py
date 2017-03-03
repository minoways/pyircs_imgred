#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
import tempfile 
from misc import *
from optparse import OptionParser

def mkcombmask(combimg, inlist, sspref='sdfr', outpref='mc', minpix=250, hsigma=3, lsigma=10, conv='block 3 3', bpm='none', sigmap='none'):
    
    # load IRAF package
    iraf.nproto()

    # chceck combined image
    if not os.access(combimg, os.R_OK):
        print >> sys.stderr, 'cannot open combined image (%s)' % combimg
        return 1

    # open input list and check if it exists
    inimg_arr = check_input(inlist, sspref)
    if isinstance(inimg_arr, int):
        return 1

    # check output images 
    outimg_arr = check_outpref(outpref, inimg_arr)
    if isinstance(outimg_arr, int):
        return 1 

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

    # check sigma map 
    if sigmap.lower() == 'none':
        tmp_sigma = ''
        tmp_mask = ''
    else:
        if os.access(sigmap, os.R_OK) == False:
            print >> sys.stderr, 'cannot read sigma map image (%s)' % sigmap
            return 1
        else:
            tmp_mask = tmp_prefix+'mask'+os.path.basename(sigmap)
            tmp_sigma = tmp_prefix+'sigma'+os.path.basename(sigmap)
            iraf.imcopy(sigmap, tmp_sigma, verbose='no')
            ret = iraf.imstat(sigmap,nclip=50,field='mode',format='no',Stdout=1)
            sigma_mode=float(ret[0])
            iraf.imreplace(tmp_sigma,sigma_mode,upper=0)
            iraf.imarith(sigmap,'*',-1,tmp_mask, verbose='no')
            iraf.imreplace(tmp_mask,1,lower=0)
            iraf.imreplace(tmp_mask,0,upper=0)
            iraf.imarith(tmp_mask,'*',-1,tmp_mask)
            iraf.imarith(tmp_mask,'+',1,tmp_mask)
    # 
    tmp_precombmask = tmp_prefix+'premask'+os.path.basename(combimg)
    tmp_combmask = tmp_prefix+'mask'+os.path.basename(combimg)
    iraf.unlearn('objmasks')
    try:
        iraf.objmasks(combimg, tmp_precombmask, omtype='boolean', masks='', convolve=conv, lsigma=lsigma, hsigma=hsigma, hdetect='yes', ldetect='yes', minpix=minpix, sigmas=tmp_sigma)
        iraf.objmasks(combimg, tmp_combmask, omtype='boolean', masks=tmp_precombmask+'[pl]', convolve=conv, lsigma=lsigma, hsigma=hsigma, hdetect='yes', ldetect='yes', minpix=minpix, sigmas=tmp_sigma)
        os.remove(tmp_precombmask)
    except:
        print >> sys.stderr, 'failed to execute iraf objmasks task for %s' % combimg
        remove_temp_all(tmp_prefix)
        return 1


    # remove the area with no sigma value 
    tmp_combmask2 = tmp_prefix+'mask2'+os.path.basename(combimg)
    if sigmap.lower() == 'none':
        iraf.imcopy(tmp_combmask+'[pl]', tmp_combmask2, verbose='no')
    else:
        iraf.imarith(tmp_combmask+'[pl]','*',tmp_mask, tmp_combmask2, verbose='no')
    os.remove(tmp_combmask)
    
    # check if input images exist and fix bad pixel if requested
    iraf.unlearn('imcopy')
    iraf.unlearn('imshift')
    for i in range(len(inimg_arr)):
        
        # get offset values from header
        im = pyfits.open(inimg_arr[i])
        try:
            dx = float(im[0].header['dx'])
            dy = float(im[0].header['dy'])
            dxc = float(im[0].header['dxc'])
            dyc = float(im[0].header['dyc'])
            nx = int(im[0].header['NAXIS1']) 
            ny = int(im[0].header['NAXIS2']) 
        except:
            print >> sys.stderr, 'failed to get offset values from the header of %s' % inimg_arr[i]
            remove_temp_all(tmp_prefix)
            return 1
        im.close()

        # shift combined mask
        tmp_shiftimg = tmp_prefix+os.path.basename(outimg_arr[i])
        iraf.imshift(tmp_combmask2, tmp_shiftimg, dxc, dyc, boundary_typ='constant', constant=0.0)
        tmp_shiftimg_reg = '%s[1:%d,1:%d]' % (tmp_shiftimg, nx, ny)
        #iraf.imcopy(tmp_shiftimg_reg, outimg_arr[i]+'[pl]',verbose='no')
        iraf.imcopy(tmp_shiftimg_reg, outimg_arr[i],verbose='no')
        if bpmex == 1:
            iraf.imarith(outimg_arr[i], '+', bpm, outimg_arr[i], verbose='no') 
        convert_maskfits_int(outimg_arr[i], outimg_arr[i])
        os.remove(tmp_shiftimg)
    
    # remove all temporary files
    remove_temp_all(tmp_prefix)
            
    return 0

if __name__=="__main__":

    usage = "usage: %prog comb_image input_list [options]"

    parser = OptionParser(usage)
    parser.add_option("--sspref", dest="sspref", type="string", default="sdfr",
                      help="prefix for sky subtracted frame (default=sdfr)")
    parser.add_option("--outpref", dest="outpref", type="string", default="mc",
                      help="prefix for output frame or output list with @ mark (default=mc)")
    parser.add_option("--minpix", dest="minpix", type="int", default=250,
                      help="minimum number of pixels in detected objects (default=250)")
    parser.add_option("--hsigma", dest="hsigma", type="float", default=3.0,
                      help="sigma threshold above sky (default=3)")
    parser.add_option("--lsigma", dest="lsigma", type="float", default=10.0,
                      help="sigma threshold below sky (default=10)")
    parser.add_option("--conv", dest="conv", type="string", default="block 3 3",
                      help="convolution kernel for objmask (default=\"block 3 3\")")
    parser.add_option("--bpm", dest="bpm", type="string", default="none",
                      help="bad pixel mask frame (default=none)")
    parser.add_option("--sigmap", dest="sigmap", type="string", default="none",
                      help="sigma map file name (default=none)")

    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()

    mkcombmask(args[0], args[1], sspref=options.sspref, outpref=options.outpref, minpix=options.minpix, hsigma=options.hsigma, lsigma=options.lsigma, conv=options.conv, bpm=options.bpm, sigmap=options.sigmap)
    
