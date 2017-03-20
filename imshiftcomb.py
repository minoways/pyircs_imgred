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

def imshiftcomb(inlist, outimg, fitgeom='shift', inpref='', objmask='none', combine='average', reject='none', fscale=False, fbase=100, fhead='F1', second=False, first_pref='sdfr', second_pref='sdf2r', indep=False, sigmap='none', expmap='none', whtmap='none', gain=5.6, ffpref=''):
    
    # check output image
    if os.access(outimg, os.R_OK):
        print >> sys.stderr, 'operation would overwrite existing image (%s)' % outimg
        return 1

    # check input image list 
    inimg_arr = check_input(inlist, inpref)
    if isinstance(inimg_arr,int):
        return 1

    # check input image list 
    if sigmap != 'none':
        ffimg_arr = check_input2(inlist, ffpref)
        if isinstance(ffimg_arr,int):
            return 1

    # check optional output image
    if sigmap != 'none':
        if os.access(sigmap, os.R_OK):
            print >> sys.stderr, 'operation would overwrite existing image (%s)' % sigmap
            return 1
    if expmap != 'none':
        if os.access(expmap, os.R_OK):
            print >> sys.stderr, 'operation would overwrite existing image (%s)' % expmap
            return 1

    if whtmap != 'none':
        if os.access(whtmap, os.R_OK):
            print >> sys.stderr, 'operation would overwrite existing image (%s)' % whtmap
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

        if os.access(dbs, os.R_OK):
            print >> sys.stderr, 'database file (%s) is already exist' % (dbs)
            return 1

        if os.access(gmp2, os.R_OK):
            print >> sys.stderr, 'modified geomap file (%s) is already exist' % (gmp2)
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

    #print 'debug'
    #print dx, max(dx), min(dx)
    #print dy, max(dy), min(dy)
        
    # check object mask
    if objmask.lower() == 'none':
        objmask = ''
    else:
        objmask_arr = check_inpref(objmask, inimg_arr)
        if isinstance(objmask_arr, int):
            return 1
                
    # independent run flag
    if indep:
        second = True
    
    # prepare for temporary file 
    tmp = tempfile.NamedTemporaryFile(suffix='', prefix='', dir='/tmp')
    tmp_prefix = tmp.name
    tmp.close()

    # calculate image median for zero shift
    iraf.unlearn('imstat')
    iraf.unlearn('mimstat')
    bgmed = []
    for i in range(len(inimg_arr)):
        if objmask == '':
            ret = iraf.imstat(inimg_arr[i], format='no', fields='midpt', nclip=50, lsigma=3., usigma=3., Stdout=1)
        else:    
            if not second:
                ret = iraf.mimstat(inimg_arr[i], imasks=objmask_arr[i]+'[pl]', format='no', fields='midpt', nclip=50, lsigma=3., usigma=3., Stdout=1)
            else:
                ret = iraf.mimstat(inimg_arr[i], imasks=objmask_arr[i], format='no', fields='midpt', nclip=50, lsigma=3., usigma=3., Stdout=1)
        if len(ret) == 1:
            bgmed.append(-1.0 * float(ret[0]))
        else:
            fout.close()
            remove_temp_all(tmp_prefix)
            print >> sys.stderr, 'failed to calculate median of the background in %s' % inimg_arr[i]
            return 1

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
    
    #print 'debug'
    #print x_size, y_size, xcmin, xcmax, ycmin, ycmax
    #print xmin, xmax, ymin, ymax

    # copy image to larger format and shift image #
    iraf.unlearn('geomap')
    iraf.unlearn('geotran')
    
    obj_list = tmp_prefix+'_obj.lst'
    if os.access(obj_list,os.R_OK):
        os.remove(obj_list)
    fobj = open(obj_list,'w')

    # for exposure time weight 
    expweight = tmp_prefix+'_exp.lst'
    if os.access(expweight,os.R_OK):
        os.remove(expweight)
    fexp = open(expweight,'w')
 
    # for zero offset 
    zeroshift = tmp_prefix+'_zeroshift.dat'
    if os.access(zeroshift, os.R_OK):
        os.remove(zeroshift)
    fzero = open(zeroshift,'w')

    # save the original fit geometry 
    fitgeom_org = fitgeom

    # preparing for the sigma list and mask
    if sigmap == 'none':
        tmp_rejmask = ''
    else:
        tmp_rejmask = tmp_prefix+'rejmask.fits'
        if os.access(tmp_rejmask, os.R_OK):
            os.remove(tmp_rejmask)
        inverse_var_list = tmp_prefix+'_var.lst'
        if os.access(inverse_var_list,os.R_OK):
            os.remove(inverse_var_list)
        finverse_var = open(inverse_var_list,'w')

    for i in range(len(inimg_arr)):

        # restore the original fit geometry 
        fitgeom = fitgeom_org 

        # geometry transformation 
        fgmp = open(gmp_arr[i])
        fgmp2 = open(gmp2_arr[i], 'w')
        nobj = 0
        for line in fgmp:
            if not line.startswith('#'):
                param = line[:-1].split()
                xref = float(param[0]) + xmin - 1
                yref = float(param[1]) + ymin - 1
                xin = float(param[2]) + xmin - 1
                yin = float(param[3]) + ymin - 1
                fgmp2.write('%.3f %.3f %.3f %.3f\n' % (xref, yref, xin, yin))
                nobj += 1
        fgmp.close()
        fgmp2.close()

        # check number of objects 
        if i == 0 and nobj == 1 and fitgeom == 'rotate':
            print 'Warning: Number of reference objects is not enought to measure the rotation' 
            print 'Warning: Only shift applied for all images'
            fitgeom = 'shift'
            fitgeom_org = 'shift' 

        if nobj == 1 and fitgeom == 'rotate':
            print 'Warning: Number of objects in %s is not enought to measure the rotation' % (inimg_arr[i])
            print 'Warning: Only shift applied for %s' % (inimg_arr[i])
            fitgeom = 'shift'

        # mapping geometry 
        iraf.geomap(gmp2_arr[i], dbs_arr[i], 1, x_size, 1, y_size, fitgeom=fitgeom, interac='no')
            
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
        img = pyfits.open(inimg_arr[i])
        if sigmap != 'none':
            ffimg = pyfits.open(ffimg_arr[i])

        # for exposure time weight
        try :
            t = float(img[0].header['EXP1TIME'])
            coadd = float(img[0].header['COADDS'])
            expt = t * coadd
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

            flux_scale = fbase / flux  
            weight = expt * (1.0 / flux_scale)**2 

        else:
            flux_scale = 1.0
            weight = expt 

        fzero.write('%f\n' % (bgmed[i]*flux_scale))
        fexp.write('%f\n' % weight)

        # object frame
        obj = np.zeros((y_size, x_size))
        obj[ymin-1:ymax,xmin-1:xmax] = img[0].data * flux_scale
        hdu = pyfits.PrimaryHDU(obj)
        obj_img = pyfits.HDUList([hdu])
        obj_img[0].header = img[0].header
        obj_img[0].header['bpm'] = msktr_fits
        obj_img[0].header['expmap'] = expt
        obj_fits = tmp_prefix+'obj'+os.path.basename(inimg_arr[i])
        objtr_fits = tmp_prefix+'objtr'+os.path.basename(inimg_arr[i])
        if os.access(obj_fits,os.R_OK):
            os.remove(obj_fits)
        obj_img.writeto(obj_fits)
        obj_img.close()
        iraf.geotran(obj_fits, objtr_fits, dbs_arr[i], gmp2_arr[i], geometr='linear', boundar='constant', constant=0)
        fobj.write('%s\n' % objtr_fits)
        img.close()

        if sigmap != 'none':
            inverse_var = np.zeros((y_size, x_size))
            inverse_var[ymin-1:ymax,xmin-1:xmax] = (np.sqrt(ffimg[0].data / (gain * expt)) * flux_scale)**-2
            hdu_var = pyfits.PrimaryHDU(inverse_var)
            inverse_var_img = pyfits.HDUList([hdu_var])
            inverse_var_img[0].header = img[0].header
            inverse_var_img[0].header['bpm2'] = '%s[*,*,%d]' % (tmp_rejmask, i+1)
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
    fobj.close()
    fexp.close()
    fzero.close()
    if sigmap != 'none':
        finverse_var.close()
    
    # combine image 
    comb_img = tmp_prefix + '_comb.fits'
    if os.access(comb_img, os.R_OK):
        os.remove(comb_img)
    if expmap == 'none':
        tmp_expmap = ''
    else:
        tmp_expmap = tmp_prefix+'expmap.fits'
        if os.access(tmp_expmap, os.R_OK):
            os.remove(tmp_expmap)

    iraf.unlearn('imcombine')
    try:
        iraf.imcombine('@'+obj_list, comb_img, sigma='', rejmask=tmp_rejmask, expmasks=tmp_expmap, combine=combine, reject=reject, masktype='!BPM', maskvalue=0.0, zero='@'+zeroshift, weight='@'+expweight, expname='EXPMAP')
    except:
        if os.access(comb_img, os.R_OK):
            os.remove(comb_img)
        if expmap != 'none':
            if os.access(tmp_expmap, os.R_OK):
                os.remove(tmp_expmap)
        iraf.imcombine('@'+obj_list, comb_img, sigma='', rejmask='', expmasks=tmp_expmap, combine=combine, reject=reject, masktype='!BPM', maskvalue=0.0, zero='@'+zeroshift, weight='@'+expweight, expname='EXPMAP')
    
    if sigmap != 'none':
        tmp_inverse_var_sum = tmp_prefix+'inverse_var.fits'
        if os.access(tmp_inverse_var_sum, os.R_OK):
            os.remove(tmp_inverse_var_sum)
        iraf.imcombine('@'+inverse_var_list, tmp_inverse_var_sum, combine='sum', reject='none', masktype='!BPM', maskvalue=0.0)
        iraf.stsdas()
        tmp_sigma = tmp_prefix+'sigma.fits'
        if os.access(tmp_sigma, os.R_OK):
            os.remove(tmp_sigma)
        iraf.imcalc(tmp_inverse_var_sum, tmp_sigma, 'sqrt(1.0/im1)', pixtype='double')
    
    # cut image
    iraf.unlearn('imcopy')
    cut_img = '%s[%d:%d,%d:%d]' % (comb_img, xcmin, xcmax, ycmin, ycmax)
    iraf.imcopy(cut_img, outimg)
    if expmap != 'none':
        cut_exp = '%s[%d:%d,%d:%d]' % (tmp_expmap, xcmin, xcmax, ycmin, ycmax)
        iraf.imcopy(cut_exp, expmap)
    if sigmap != 'none':
        cut_sigma = '%s[%d:%d,%d:%d]' % (tmp_sigma, xcmin, xcmax, ycmin, ycmax)
        iraf.imcopy(cut_sigma, sigmap)
        if whtmap != 'none':
            cut_wht = '%s[%d:%d,%d:%d]' % (tmp_inverse_var_sum, xcmin, xcmax, ycmin, ycmax)
            iraf.imcopy(cut_wht, whtmap)
    
    # delete temporary object files
    remove_temp_all(tmp_prefix+'obj')
    os.remove(obj_list)
    os.remove(comb_img)

    # record relative offset between input images and combined image and rotation
    for i in range(len(inimg_arr)):
        im = pyfits.open(inimg_arr[i],mode='update')

        if second:
            
            if indep:

                # calculate offset
                dxc = xcmin - xmin - dx[i]
                dyc = ycmin - ymin - dy[i]

                # retrieve rotation
                rot = 0.0
                fdbs = open(dbs_arr[i])
                for line in fdbs:
                    param = line[:-1].split()
                    if param[0] == 'xrotation':
                        rot = float(param[1])

                if rot > 180.0:
                    rot = rot - 360.0

                im[0].header['dx'] = dx[i]
                im[0].header['dy'] = dy[i]
                im[0].header['dxc'] = dxc
                im[0].header['dyc'] = dyc
                im[0].header['rotation'] = rot
                
            else:
                # check number of objects in the geomap file
                nobj = 0
                fgmp = open(gmp_arr[0])
                for line in fgmp:
                    nobj += 1
            
                im1 = pyfits.open(inimg_arr[i].replace(second_pref, first_pref),mode='update')
                for j in range(nobj):
                    key = 'XC%d' % (j+1)
                    im[0].header[key] = float(im1[0].header[key])
                    key = 'YC%d' % (j+1)
                    im[0].header[key] = float(im1[0].header[key])
                    key = 'PEAK%d' % (j+1)
                    im[0].header[key] = float(im1[0].header[key])
                    key = 'FWHM%d' % (j+1)
                    im[0].header[key] = float(im1[0].header[key])
                key = 'DX'
                im[0].header[key] = float(im1[0].header[key])
                key = 'DY'
                im[0].header[key] = float(im1[0].header[key])
                key = 'DXC'
                im[0].header[key] = float(im1[0].header[key])
                key = 'DYC'
                im[0].header[key] = float(im1[0].header[key])
                key = 'ROTATION'
                im[0].header[key] = float(im1[0].header[key])
                im1.close()    
            
        else:

            # calculate offset
            dxc = xcmin - xmin - dx[i]
            dyc = ycmin - ymin - dy[i]

            # retrieve rotation
            rot = 0.0
            fdbs = open(dbs_arr[i])
            for line in fdbs:
                param = line[:-1].split()
                if param[0] == 'xrotation':
                    rot = float(param[1])

            if rot > 180.0:
                rot = rot - 360.0

            im[0].header['dx'] = dx[i]
            im[0].header['dy'] = dy[i]
            im[0].header['dxc'] = dxc
            im[0].header['dyc'] = dyc
            im[0].header['rotation'] = rot

        im.close()

    # remove all temporary files
    remove_temp_all(tmp_prefix)

    return 0 

if __name__=="__main__":

    usage = "usage: %prog input_list output_image [options]"

    parser = OptionParser(usage)
    parser.add_option("--fitgeom", dest="fitgeom", type="choice", default="shift",
                      choices=["shift", "rotate"],
                      help="Fitting geometry (shift|rotate, default=shift)")
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=none)")
    parser.add_option("--objmask", dest="objmask", type="string", default="none",
                      help="prefix for mask image, image list with \'@\', or header keyword with \':\'")
    parser.add_option("--combine", dest="combine", type="choice", default="average",
                      choices=["average","median","sum"],
                      help="type of combine operation (average|median|sum)")
    parser.add_option("--reject", dest="reject", type="choice", default="sigclip",
                      choices=["none", "minmax", "ccdclip", "crreject", "sigclip", "avsigclip", "pclipor"],
                      help="type of rejection (none|minmax|ccdclip|crreject|sigclip|avsigclip|pclipor, default=sigclip)")
    parser.add_option("--fscale", dest="fscale", action="store_true", default=False,
                      help="Scaling with the reference object flux (default=False)")    
    parser.add_option("--fbase", dest="fbase", type="float", default=100.0,
                      help="Baseline flux for the flux scaling (default=100.0)")    
    parser.add_option("--fhead", dest="fhead", type="string", default="F1",
                      help="Header keyword for the reference object flux (default=F1)")    
    parser.add_option("--second", dest="second", action="store_true", default=False,
                      help="Is this second cycle subtraction? (default=False)")    
    parser.add_option("--first_pref", dest="first_pref", type="string", default='sdfr',
                      help="Prefix of the first sky subtracted images (default=sdfr)")    
    parser.add_option("--second_pref", dest="second_pref", type="string", default='sdf2r',
                      help="Prefix of the second sky subtracted images (default=sdf2r)")    
    parser.add_option("--indep", dest="indep", action="store_true", default=False,
                      help="Flag for independent command execution (default=False)")
    parser.add_option("--gain", dest="gain", type="float", default=5.6,
                      help="Detector gain [e/ADU] (default=5.6)")    
    parser.add_option("--sigmap", dest="sigmap", type="string", default="none",
                      help="Sigma map file name (default=none)")
    parser.add_option("--expmap", dest="expmap", type="string", default="none",
                      help="Exposure map file name (default=none)")
    parser.add_option("--whtmap", dest="whtmap", type="string", default="none",
                      help="Weight map file name (default=none)")
    parser.add_option("--ffpref", dest="ffpref", type="string", default="",
                      help="prefix for flat frame (to be replaced from the original prefix in the input list) (default=none)")

    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        sys.exit()

    imshiftcomb(args[0], args[1], fitgeom=options.fitgeom, inpref=options.inpref, objmask=options.objmask, combine=options.combine, reject=options.reject, fscale=options.fscale, fbase=options.fbase, fhead=options.fhead, second=options.second, first_pref=options.first_pref, second_pref=options.second_pref, indep=options.indep, expmap=options.expmap, sigmap=options.sigmap, whtmap=options.whtmap, gain=options.gain, ffpref=options.ffpref)
    
