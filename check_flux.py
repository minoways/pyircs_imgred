#!/usr/bin/env python

import os,sys
import pyfits
from misc import *
from pyraf import iraf
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 
from scipy.ndimage import median_filter

def check_flux(inlist, outdata, inpref='', rap=50, salgor='median', annu=50, dannu=5, skyconst=0.0, xhead='XC1', yhead='YC1', fhead='F1'):

    # check input image list 
    inimg_arr = check_input(inlist, inpref)
    if isinstance(inimg_arr,int):
        return 1

    # check output
    if os.access(outdata, os.R_OK):
        print >> sys.stderr, 'operation would overwrite existing data (%s)' % outdata
        return 1
    else:
        fout = open(outdata,'w')
    
    # unlearn IRAF commands 
    iraf.unlearn('hedit')
    
    frame_id = []
    flux_arr = []
    for i in range(len(inimg_arr)):

        # read data into memory
        im = pyfits.open(inimg_arr[i], mode='update')
        im_data = im[0].data
        xc = float(im[0].header[xhead])
        yc = float(im[0].header[yhead])

        x = np.ones(im_data.shape[0]*im_data.shape[1]).reshape(im_data.shape[0],im_data.shape[1])* np.arange(im_data.shape[1]) + 1
        y = (np.ones(im_data.shape[1]*im_data.shape[0]).reshape(im_data.shape[1],im_data.shape[0]) * np.arange(im_data.shape[0]) + 1).T
        r = np.sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc))

        # define the area to measure the flux
        phot_idx = np.where(r<=rap)
        th_idx = np.where(r<=annu+dannu)

        # remove bad pixels 
        blurred = median_filter(im_data, size=3)
        diff = im_data - blurred
        threshold = 3.0 * np.nanstd(diff[th_idx])

        # find the hot pixels 
        hotpix = np.nonzero((np.abs(diff[1:-1,1:-1])>threshold))
        hotpix = np.array(hotpix) + 1

        # fix the hot pixels 
        fixed_im_data = np.copy(im_data)
        for y,x in zip(hotpix[0], hotpix[1]):
            fixed_im_data[y,x] = blurred[y,x]

        # sky estimation
        if salgor == "constant":
            sky = skyconst
        else:
            sky_idx = np.logical_and(r>=annu, r<annu+dannu)
            sky = np.median(fixed_im_data[sky_idx])

        # calculate flux 
        flux = np.sum(fixed_im_data[phot_idx]) - len(fixed_im_data[phot_idx]) * sky 
        
        # record flux 
        im[0].header[fhead] = flux
        fout.write('%s %d %f\n' % (inimg_arr[i], i, flux))
        frame_id.append(i)
        flux_arr.append(flux)

        # close image handler 
        im.close()

    # close file handler
    fout.close()

    # plot result
    plt.ion()
    plt.figure()
    plt.plot(frame_id, flux_arr,'*')
    plt.xlabel('Frame number')
    plt.ylabel('Flux [ADU/sec]')
    plt.grid()
    plt.show()
    _ = raw_input("Hit enter to finish the script:")
    
    return 0
    

if __name__=="__main__":

    usage = "usage: %prog input_list output_data [options]"

    parser = OptionParser(usage)
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=none)")
    parser.add_option("--rap", dest="rap", type="int", default=50,
                      help="Aperture radius in pix (default=50)")
    parser.add_option("--salgor", dest="salgor", type="choice", default="median",
                      choices=["median","constant"],
                      help="sky fitting algroithm (median|constant)")
    parser.add_option("--annu", dest="annu", type="int", default=50,
                      help="Inner radius of sky annulus in pix (default=50)")
    parser.add_option("--dannu", dest="dannu", type="int", default=5,
                      help="Width of sky annulus in pix (default=5)")
    parser.add_option("--skyconst", dest="skyconst", type="float", default=0.0,
                      help="Constant sky value per pixel (default=0.0)")
    parser.add_option("--xhead", dest="xhead", type="string", default="XC1",
                      help="Header keyword for the X coordinate of the object (default=XC1)")
    parser.add_option("--yhead", dest="yhead", type="string", default="YC1",
                      help="Header keyword for the Y coordinate of the object (default=YC1)")
    parser.add_option("--fhead", dest="fhead", type="string", default="F1",
                      help="Header keyword for recording the measured flux (default=F1)")
    
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()

    ret = check_flux(args[0],args[1], inpref=options.inpref, rap=options.rap, salgor=options.salgor, annu=options.annu, dannu=options.dannu, skyconst=options.skyconst, xhead=options.xhead, yhead=options.yhead, fhead=options.fhead)

 
    
    
