#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
import tempfile 
from zcombine import * 
from misc import *
import math 
import numpy as np
from optparse import OptionParser
import shutil 

def get_large_region(nx, ny, dx, dy):

    xpad = []
    dx_min = min(dx)
    if dx_min > 0:
        xpad.append(0.0)
    else :
        xpad.append(math.fabs(dx_min))
    dx_max = max(dx)
    if dx_max < 0:
        xpad.append(0.0)
    else:
        xpad.append(dx_max)
    xpad_size = int(max(xpad)) + 1

    ypad = []
    dy_min = min(dy)
    if dy_min > 0:
        ypad.append(0.0)
    else:
        ypad.append(math.fabs(dy_min))
    dy_max = max(dy)
    if dy_max < 0:
        ypad.append(0.0)
    else:
        ypad.append(dy_max)
    ypad_size = int(max(ypad))+1

    x_size = nx + xpad_size + xpad_size
    y_size = ny + ypad_size + ypad_size
    
    xcmin = int((x_size-nx)/2 + 1 + dx_min)
    xcmax = int((x_size-nx)/2 + nx + dx_max)
    ycmin = int((y_size-ny)/2 + 1 + dy_min)
    ycmax = int((y_size-ny)/2 + ny + dy_max)

    return x_size, y_size, xcmin, xcmax, ycmin, ycmax

def sigmap(inlist, sigmap, expmap='none', whtmap='none', inpref='', ffpref='', objmask='none', reject='sigclip', fscale=False, fbase=100, fhead='F1', gain=5.6):
    
    # check output image
    if os.access(sigmap, os.R_OK):
        print >> sys.stderr, 'operation would overwrite existing image (%s)' % sigmap
        return 1

    if os.access(expmap, os.R_OK):
        print >> sys.stderr, 'operation would overwrite existing image (%s)' % expmap
        return 1

    if os.access(whtmap, os.R_OK):
        print >> sys.stderr, 'operation would overwrite existing image (%s)' % whtmap
        return 1

    # check input image list 
    inimg_arr = check_input(inlist, inpref)
    if isinstance(inimg_arr,int):
        return 1

    # check input image list 
    ffimg_arr = check_input2(inlist, ffpref)
    if isinstance(ffimg_arr,int):
        return 1
    
    # get array size 
    im = pyfits.open(inimg_arr[0])
    nx = im[0].header['NAXIS1'] 
    ny = im[0].header['NAXIS2'] 
    im.close()
    
    # check geomap data
    dx = []
    dy = []
    gmp_arr = []
    gmp2_arr = []
    dbs_arr = []
    for i in range(len(inimg_arr)):
        fname,ext = os.path.splitext(inimg_arr[i])
        gmp = fname + '.gmp'
        gmp2 = fname + '.gmp2'
        dbs = fname + '.dbs'

        gmp_arr.append(gmp)
        gmp2_arr.append(gmp2)
        dbs_arr.append(dbs)

        if not os.access(gmp, os.R_OK):
            print >> sys.stderr, 'geomap file (%s) does not exist' % (gmp)
            return 1

        if not os.access(dbs, os.R_OK):
            print >> sys.stderr, 'database file (%s) does not exist' % (dbs)
            return 1

        if not os.access(gmp2, os.R_OK):
            print >> sys.stderr, 'modified geomap file (%s) does not exist' % (gmp2)
            return 1

        fgmp = open(gmp)
        nl = 1
        dx_ave = 0.0
        dy_ave = 0.0
        for line in fgmp:
            if not line.startswith('#'):
                param = line[:-1].split()
                if len(param) != 4:
                    print >> sys.stderr, 'Invalid format in line %d of %s: %s' % (nl, gmp, line[:-1])
                    fgmp.close()
                    return 1
                else:
                    if isfloat(param[0]) == False or isfloat(param[1]) == False or isfloat(param[2]) == False or isfloat(param[3]) == False:
                        print >> sys.stderr, 'failed to decode line %d of %s: %s' % (nl, gmp, line[:-1])
                        fgmp.close()
                        return 1
                    else :
                        dx_ave += float(param[0]) - float(param[2])
                        dy_ave += float(param[1]) - float(param[3])

                nl += 1
        #print inimg_arr[i],nl
        dx.append(dx_ave / (nl-1))
        dy.append(dy_ave / (nl-1))

    if len(inimg_arr) != len(dx):
        print >> sys.stderr, 'number of input images does not match with that of offsets'
        return 1
       
    # check object mask
    if objmask.lower() == 'none':
        objmask = ''
    else:
        objmask_arr = check_inpref(objmask, inimg_arr)
        if isinstance(objmask_arr, int):
            return 1

    # prepare for temporary file 
    tmp = tempfile.NamedTemporaryFile(suffix='', prefix='', dir='/tmp')
    tmp_prefix = tmp.name
    tmp.close()

    # get large array size and combined image size 
    ret = get_large_region(nx, ny, dx, dy)
    if len(ret) != 6:
        print >> sys.stderr, 'failed to get large array size' 
        return 1
    
    x_size = ret[0]
    y_size = ret[1]
    xcmin = ret[2]
    xcmax = ret[3]
    ycmin = ret[4]
    ycmax = ret[5]
    
    # calculate image region in the large format
    xmin = int((x_size - nx)/2)+1
    xmax = nx + int((x_size - nx)/2)
    ymin = int((y_size - ny)/2)+1
    ymax = ny + int((y_size - ny)/2)
    
    # copy image to larger format and shift image #
    iraf.unlearn('geomap')
    iraf.unlearn('geotran')
    
    # for exposure time weight and flux scaling
    expt_arr = []
    flux_scale_arr = []
    for i in range(len(inimg_arr)):

        # load original frame
        img = pyfits.open(inimg_arr[i])
        
        # for exposure time weight
        try :
            t = float(img[0].header['EXP1TIME'])
            coadd = float(img[0].header['COADDS'])
            expt_arr.append(t * coadd)
        except KeyError:
            print >> sys.stderr, 'can not read exposure time from the header of %s' % inimg_arr[i]
            img.close()
            return 1

        # for flux scaling and weight 
        if fscale:
            try: 
                flux = float(img[0].header[fhead])
            except KeyError:
                print >> sys.stderr, 'can not read flux keyword (%s) from the header of %s' % (fhead, inimg_arr[i])
                img.close()
                return 1
            flux_scale_arr.append(fbase / flux)  
        else:
            flux_scale_arr.append(1.0)
        img.close()
    
    # preparing weighted variance map for each image
    inverse_var_list = tmp_prefix+'_obj.lst'
    if os.access(inverse_var_list,os.R_OK):
        os.remove(inverse_var_list)
    finverse_var = open(inverse_var_list,'w')

    for i in range(len(inimg_arr)):

        # mask frame
        msk = np.ones((y_size,x_size))
        msk[ymin-1:ymax,xmin-1:xmax] = 0
        hdu = pyfits.PrimaryHDU(msk)
        msk_img = pyfits.HDUList([hdu])
        msk_fits = tmp_prefix+'mask'+os.path.basename(inimg_arr[i])
        msktr_fits = tmp_prefix+'masktr'+os.path.basename(inimg_arr[i])
        if os.access(msk_fits,os.R_OK):
            os.remove(msk_fits)
        if os.access(msktr_fits,os.R_OK):
            os.remove(msktr_fits)
        msk_img.writeto(msk_fits)
        msk_img.close()
        
        # transform mask geometry
        iraf.geotran(msk_fits, msktr_fits, dbs_arr[i], gmp2_arr[i], geometr='linear', boundar='constant', constant=1)
        os.remove(msk_fits)
        convert_maskfits_int(msktr_fits, msktr_fits)
        
        # load original frame
        ffimg = pyfits.open(ffimg_arr[i])
        img = pyfits.open(inimg_arr[i])
        # object frame
        inverse_var = np.zeros((y_size, x_size))
        #print np.median(ffimg[0].data), ffimg_arr[i], gain, expt_arr[i], flux_scale_arr[i]
        #print np.median(np.sqrt(ffimg[0].data / (gain * expt_arr[i]))), flux_scale_arr[i]
        inverse_var[ymin-1:ymax,xmin-1:xmax] = (np.sqrt(ffimg[0].data / (gain * expt_arr[i])) * flux_scale_arr[i])**-2
        hdu = pyfits.PrimaryHDU(inverse_var)
        inverse_var_img = pyfits.HDUList([hdu])
        inverse_var_img[0].header = img[0].header
        inverse_var_img[0].header['bpm'] = msktr_fits
        inverse_var_img[0].header.update('EXPTIME', expt_arr[i])
        #inverse_var_img[0].header.update('MASKSCAL', expt_arr[i])
        #inverse_var_img[0].header.update('MASKZERO', expt_arr[i])
        #print 'EXPT = %f' % (expt_arr[i])
        inverse_var_fits = tmp_prefix+'var'+os.path.basename(inimg_arr[i])
        inverse_vartr_fits = tmp_prefix+'vartr'+os.path.basename(inimg_arr[i])
        if os.access(inverse_var_fits,os.R_OK):
            os.remove(inverse_var_fits)
        if os.access(inverse_vartr_fits,os.R_OK):
            os.remove(inverse_vartr_fits)
        inverse_var_img.writeto(inverse_var_fits)
        inverse_var_img.close()
        iraf.geotran(inverse_var_fits, inverse_vartr_fits, dbs_arr[i], gmp2_arr[i], geometr='linear', boundar='constant', constant=0)
        finverse_var.write('%s\n' % inverse_vartr_fits)
        ffimg.close()

    # close file handlers
    finverse_var.close()

    # sum weighted variance images 
    tmp_inverse_var_sum = tmp_prefix+'inverse_var.fits'
    if os.access(tmp_inverse_var_sum, os.R_OK):
        os.remove(tmp_inverse_var_sum)
    tmp_sigma = tmp_prefix+'sigma.fits'
    if os.access(tmp_sigma, os.R_OK):
        os.remove(tmp_sigma)
    tmp_exp = tmp_prefix+'exp.fits'
    if os.access(tmp_exp, os.R_OK):
        os.remove(tmp_exp)
    
    if expmap != 'none':
        #iraf.hselect('@'+inverse_var_list,"$I,EXPTIME,MASKSCAL,MASKZERO","yes")
        iraf.imcombine('@'+inverse_var_list, tmp_inverse_var_sum, expmasks=tmp_exp, combine='sum', reject=reject, masktype='!BPM', maskvalue=0.0, expname='EXPTIME')
    else:
        iraf.imcombine('@'+inverse_var_list, tmp_inverse_var_sum, combine='sum', reject=reject, masktype='!BPM', maskvalue=0.0)

    # calculate sigma 
    iraf.stsdas()
    iraf.imcalc(tmp_inverse_var_sum, tmp_sigma, 'sqrt(1.0/im1)', pixtype='double')
        
    # cut image
    iraf.unlearn('imcopy')
    cut_sig = '%s[%d:%d,%d:%d]' % (tmp_sigma, xcmin, xcmax, ycmin, ycmax)
    iraf.imcopy(cut_sig, sigmap)
    if expmap != 'none':
        cut_exp = '%s[%d:%d,%d:%d]' % (tmp_exp, xcmin, xcmax, ycmin, ycmax)
        iraf.imcopy(cut_exp, expmap)

    # calc weight map
    if whtmap != 'none':
        cut_wht = '%s[%d:%d,%d:%d]' % (tmp_inverse_var_sum, xcmin, xcmax, ycmin, ycmax)
        iraf.imcopy(cut_wht, whtmap)
    
    # delete temporary object files
    os.remove(inverse_var_list)

    # remove all temporary files
    #remove_temp_all(tmp_prefix)

    return 0 

if __name__=="__main__":

    usage = "usage: %prog input_list sigma_map [options]"

    parser = OptionParser(usage)
    parser.add_option("--expmap", dest="expmap", type="string", default="none",
                      help="exposure map file name (default=none)")
    parser.add_option("--whtmap", dest="whtmap", type="string", default="none",
                      help="weight map file name (default=none)")
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=none)")
    parser.add_option("--ffpref", dest="ffpref", type="string", default="",
                      help="prefix for flat frame (to be replaced from the original prefix in the input list) (default=none)")
    parser.add_option("--objmask", dest="objmask", type="string", default="none",
                      help="prefix for mask image, image list with \'@\', or header keyword with \':\'")
    parser.add_option("--reject", dest="reject", type="choice", default="sigclip",
                      choices=["none", "minmax", "ccdclip", "crreject", "sigclip", "avsigclip", "pclipor"],
                      help="type of rejection (none|minmax|ccdclip|crreject|sigclip|avsigclip|pclipor, default=sigclip)")
    parser.add_option("--fscale", dest="fscale", action="store_true", default=False,
                      help="Scaling with the reference object flux (default=False)")    
    parser.add_option("--fbase", dest="fbase", type="float", default=100.0,
                      help="Baseline flux for the flux scaling (default=100.0)")    
    parser.add_option("--fhead", dest="fhead", type="string", default="F1",
                      help="Header keyword for the reference object flux (default=F1)")    
    parser.add_option("--gain", dest="gain", type="float", default=5.6,
                      help="Detector gain [e/ADU] (default=5.6)")    

    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        sys.exit()

    sigmap(args[0], args[1], expmap=options.expmap, whtmap=options.whtmap, inpref=options.inpref, ffpref=options.ffpref, objmask=options.objmask, reject=options.reject, fscale=options.fscale, fbase=options.fbase, fhead=options.fhead, gain=options.gain)
    
