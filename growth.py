#!/usr/bin/env python

import os,sys
import pyfits
from misc import *
from ds9 import ds9 
from pyraf import iraf
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 
from scipy import stats

def growth(inimg, outdata, rmin=1, rmax=50, rbin=1, salgor='median', annu=50, dannu=5, skyconst=0.0):

    # check input image
    if not os.access(inimg, os.R_OK):
        print >> sys.stderr, 'cannot read input image (%s)' % inimg
        return 1

    # check output
    if os.access(outdata, os.R_OK):
        print >> sys.stderr, 'operation would overwrite existing data (%s)' % outdata
        return 1
    else:
        fout = open(outdata,'w')
    
    # open ds9 and display image
    d = ds9()
    d.set('zoom to fit')

    # unlearn IRAF commands 
    iraf.unlearn('display')
    iraf.unlearn('rimexam')
    iraf.unlearn('imexam')
    iraf.unlearn('hedit')
    
    # display image
    d.set('regions delete all')
    iraf.display(inimg, 1)

    # imexam
    print '\nSelect Object --> type a , Quit --> type q on the ds9 image'
    print 'Caution: do not type any key except a or q'
    ret = iraf.imexam(inimg, 1, Stdout=1)
    
    # display result
    for j in range(len(ret)):
        print ret[j]
    print '\n'
            
    # parse results
    param1 = ret[len(ret)-2].split()
    param2 = ret[len(ret)-1].split()
    if len(param1) == 4 and len(param2) == 11:
        if isfloat(param1[0]):
            xobj = float(param1[0])
        if isfloat(param1[0]):
            yobj = float(param1[1])
        if isfloat(param2[4]):
            peak = float(param2[4])
        if isfloat(param2[9]):
            fwhm = float(param2[9])
        else :
            print >> sys.stderr, 'failed to pick up object in %s' % inimg
            return 1

    # read data into memory
    im = pyfits.open(inimg)
    im_data = im[0].data
    im.close()
    
    x = np.ones(im_data.shape[0]*im_data.shape[1]).reshape(im_data.shape[0],im_data.shape[1])* np.arange(im_data.shape[1]) + 1
    y = (np.ones(im_data.shape[1]*im_data.shape[0]).reshape(im_data.shape[1],im_data.shape[0]) * np.arange(im_data.shape[0]) + 1).T
    r = np.sqrt((x-xobj)*(x-xobj) + (y-yobj)*(y-yobj))

    # sky estimation
    if salgor == "constant":
        sky = skyconst
    else:
        sky_idx = np.logical_and(r>=annu, r<annu+dannu)
        sky = np.median(im_data[sky_idx])
    
    # growth curve
    rap_arr = []
    pix_arr = []
    flux_arr = []
    rap = rmin
    while rap <= annu + dannu:
        phot_idx = np.where(r<=rap)
        total_pix = len(im_data[phot_idx])
        flux = np.sum(im_data[phot_idx]) - total_pix * sky 
        pix_arr.append(total_pix)
        rap_arr.append(rap)
        flux_arr.append(flux)
        rap = rap + rbin
        
    # update sky level
    if salgor != "constant":
        npix1 = int(annu**2 * math.pi)
        npix2 = int((annu + dannu)**2 * math.pi)
        f_idx = np.logical_and(np.array(pix_arr)>npix1, np.array(pix_arr)<=npix2)
        aa, bb, rr, _, _ = stats.linregress(np.array(pix_arr)[f_idx], np.array(flux_arr)[f_idx])
        sky = sky + aa 

        print 'Sky level = %f [ADU/pix]' % sky

    # growth curve
    rap_arr = []
    pix_arr = []
    flux_arr = []
    rap = rmin
    while rap <= rmax:
        phot_idx = np.where(r<=rap)
        total_pix = len(im_data[phot_idx])
        flux = np.sum(im_data[phot_idx]) - total_pix * sky 
        fout.write('%d %f\n' % (rap, flux))
        pix_arr.append(total_pix)
        rap_arr.append(rap)
        flux_arr.append(flux)
        rap = rap + rbin
    fout.close()
    
    return rap_arr, flux_arr
    

if __name__=="__main__":

    usage = "usage: %prog input_image output_data [options]"

    parser = OptionParser(usage)    
    parser.add_option("--rmin", dest="rmin", type="int", default=1,
                      help="minimum radius in pix (default=1)")
    parser.add_option("--rmax", dest="rmax", type="int", default=50,
                      help="maximum radius in pix (default=50)")
    parser.add_option("--rbin", dest="rbin", type="int", default=1,
                      help="Resolution of growth curve in pix (default=1)")
    parser.add_option("--salgor", dest="salgor", type="choice", default="median",
                      choices=["median","constant"],
                      help="sky fitting algroithm (median|constant)")
    parser.add_option("--annu", dest="annu", type="int", default=50,
                      help="Inner radius of sky annulus in pix (default=50)")
    parser.add_option("--dannu", dest="dannu", type="int", default=5,
                      help="Width of sky annulus in pix (default=5)")
    parser.add_option("--skyconst", dest="skyconst", type="float", default=0.0,
                      help="Constant sky value per pixel (default=0.0)")
    
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()

    ret = growth(args[0],args[1], rmin=options.rmin, rmax=options.rmax, rbin=options.rbin, salgor=options.salgor, annu=options.annu, dannu=options.dannu, skyconst=options.skyconst)

    if not isinstance(ret, int):

        # measure total flux
        total_flux = ret[1][len(ret[1]) - 1]
    
        print 'Total flux = %.2f\n' % (total_flux)

        # plot 
        plt.ion()
        plt.figure()
        plt.plot(ret[0], ret[1],'*')
        plt.xlabel('Aperture radiux [pix]')
        plt.ylabel('Encircled flux')
        plt.grid()
        plt.show()
        _ = raw_input("Hit enter to finish the script.")        

 
    
    
