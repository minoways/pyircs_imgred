#!/usr/bin/env python

import os,sys
import pyfits
from ds9 import ds9 
from misc import *
from optparse import OptionParser
import math
import numpy as np
from pyraf import iraf
import matplotlib.pyplot as plt 
import tempfile
from zcombine import * 
from ds9 import ds9

def frcheck_auto(infile, inpref='', clipsig=3.0):
    
    # open input list and check if it exists
    inimg_arr = check_input(infile, inpref)
    if isinstance(inimg_arr, int):
        return 1

    # open list and backup original file
    f=open(infile)
    orgfile = infile+".org"
    f_org = open(orgfile,'w')
    for line in f:
        f_org.write(line)
    # close list 
    f.close()
    f_org.close()

    # prefix for temoprary images
    tmp = tempfile.NamedTemporaryFile(suffix='', prefix='', dir='/tmp')
    tmp_prefix = tmp.name
    tmp.close()

    # normalize by coadd,ndr,exptime and combine 
    nframe = 0
    tmp_raw_list = '%s_raw.lst' % tmp_prefix
    tmp_list = []
    f = open(tmp_raw_list, 'w')
    x = []
    label = []
    for i in range(len(inimg_arr)):
        tmp_fits = '%s_%s' % (tmp_prefix,os.path.basename(inimg_arr[i]))
        f.write(tmp_fits+'\n')
        im = pyfits.open(inimg_arr[i])
        expt = float(im[0].header['EXP1TIME'])
        coadd = int(im[0].header['COADDS'])
        ndr = int(im[0].header['NDR'])
        scale = expt * coadd * ndr 
        im.close()
        iraf.imarith(inimg_arr[i],'/',scale, tmp_fits, verbose='no')
        x.append(i)
        tmp_list.append(tmp_fits)
        label.append(os.path.basename(inimg_arr[i]).replace('.fits',''))
    f.close()

    tmp_comb = '%s_comb.fits' % tmp_prefix
    ret = zcombine(tmp_raw_list, tmp_comb, combine='median', reject='sigclip')
    if ret != 0:
        remove_temp_all(tmp_prefix)
        return 1

    # divided by temporary combined image
    try :
        iraf.imarith('@'+tmp_raw_list, '/', tmp_comb, '@'+tmp_raw_list, verbose='no')
    except:
        remove_temp_all(tmp_prefix)
        return 1

    # calculate standard deviation of the divided images
    try:
        iraf.unlearn('imstat')
        stddev = iraf.imstat('@'+tmp_raw_list, format='no', fields='stddev', nclip=50, Stdout=1)
    except:
        remove_temp_all(tmp_prefix)
        return 1

    # calculate mean and sigma of the stddev of the image
    y = map(float, stddev)
    mean, sigma = meanclip(np.array(y, dtype=float))

    # flag bad frames
    frame_flag = []
    for i in range(len(stddev)):
        if math.fabs(y[i] - mean) < clipsig * sigma:
            frame_flag.append(0)
        else:
            frame_flag.append(1)

    # 
    xr = []
    yr = []
    for i in range(len(frame_flag)):
        if frame_flag[i] == 1:
            xr.append(x[i])
            yr.append(y[i])

    plt.figure(num=None, figsize=(15, 10), dpi=80, edgecolor='k')
    plt.ion()
    plt.plot(x,y,'ro-')
    plt.plot(xr,yr,'x',linewidth=10,markersize=10,label='bad')
    plt.xticks(x,label, rotation='vertical')
    plt.margins(0.2)
    plt.subplots_adjust(bottom=0.20)

    # open ds9 
    d = ds9()
    d.set('frame delete all')
    d.set('tile')
    d.set('tile mode grid')
    for i in range(len(x)):
        if i in xr:
            iraf.hedit(tmp_list[i],'OBJECT',os.path.basename(inimg_arr[i])+'(NG)',delete='yes',verify='no',add='yes')
        else:
            iraf.hedit(tmp_list[i],'OBJECT',os.path.basename(inimg_arr[i])+'(OK)',delete='yes',verify='no',add='yes')

        d.set('frame %d' % (i+1))
        d.set('file %s' % tmp_list[i])
        d.set('scale zscale')
        if i in xr:
            d.set('cmap bb')
        else:
            d.set('cmap grey')

    for i in range(len(x)):
        d.set('frame %d' % (i+1))
        d.set('zoom to fit')
    
    # remove temporary file
    remove_temp_all(tmp_prefix)

    # check result
    check = raw_input('Is this result ok? (yes/no)')
    if check.lower() == '' or check.lower() == 'yes':
        check = 'ok'
    else:
        check = 'ng'

    if check == 'ng':
        return 0

    #<--- output bad removed list --->#
    # open list
    try:
        f=open(orgfile)
    except IOError:
        print >> sys.stderr, 'cannot read original file %s' % (orgfile)
    f_out = open(infile, 'w')
    ngfile = infile+".ng"
    fng_out = open(ngfile,'w')
    i = 0
    for line in f:
        if not line.startswith('#'):
            if frame_flag[i] == 0:
                f_out.write('%s\n' % inimg_arr[i])
            else :
                f_out.write('#%s\n' % inimg_arr[i])
                fng_out.write('%s\n' % inimg_arr[i])
            i += 1
        else:
            f_out.write(line)
            
    # close list 
    f.close()
    f_out.close()
    fng_out.close()

    return 0


if __name__=="__main__":
    
    usage = "usage: %prog file_list [options]"
    parser = OptionParser(usage)
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=none)")
    parser.add_option("--clipsig", dest="clipsig", type="float", default=3.0,
                      help="clipping sigma threshold (default=3.0)")

    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit()

    frcheck_auto(sys.argv[1], inpref=options.inpref, clipsig=options.clipsig)
