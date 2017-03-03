#!/usr/bin/env python

import os,sys
import pyfits
from misc import *
from ds9 import ds9 
from pyraf import iraf
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 
from scipy.optimize import curve_fit
import math 

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.0*sigma**2))

def calcnoise(inimg, rap, xmin=None, xmax=None, ymin=None, ymax=None, skywidth=5.0, maxap=500, nhist=30, hz1=None, hz2=None):

    # check input image
    if not os.access(inimg, os.R_OK):
        print >> sys.stderr, 'cannot read input image (%s)' % inimg
        return 1

    im = pyfits.open(inimg)
    im_data = im[0].data
    nx = int(im[0].header['NAXIS1']) 
    ny = int(im[0].header['NAXIS2'])  
    im.close()

    # check range
    if xmin == None:
        xmin = 1 
    if xmax == None:
        xmax = nx
    if ymin == None:
        ymin = 1
    if ymax == None:
        ymax = ny

    # check maximum number of apertures 
    xwidth = xmax - xmin + 1 
    ywidth = ymax - ymin + 1 
    buf = int(rap) * 2 
    x_range = xwidth - 2 * buf
    y_range = ywidth - 2 * buf 
    area = x_range * y_range
    apnum = int(area / (rap*rap*math.pi))
    if apnum > maxap:
        apnum = maxap
        
    # generate randome coordinates 
    xrand = np.random.rand(apnum) * x_range + (xmin + buf)
    yrand = np.random.rand(apnum) * x_range + (xmin + buf)

    # calculate flux
    flux_arr = []
    x = np.ones(im_data.shape[0]*im_data.shape[1]).reshape(im_data.shape[0],im_data.shape[1])* np.arange(im_data.shape[1]) + 1
    y = (np.ones(im_data.shape[1]*im_data.shape[0]).reshape(im_data.shape[1],im_data.shape[0]) * np.arange(im_data.shape[0]) + 1).T
    for i in range(apnum):
        x1 = int(xrand[i] - (buf + skywidth))
        x2 = int(xrand[i] + (buf + skywidth))
        y1 = int(yrand[i] - (buf + skywidth))
        y2 = int(yrand[i] + (buf + skywidth))
        im_sub = im_data[y1-1:y2,x1-1:x2]
        x_sub = x[y1-1:y2,x1-1:x2]
        y_sub = y[y1-1:y2,x1-1:x2]
        r = np.sqrt((x_sub-xrand[i])*(x_sub-xrand[i]) + (y_sub-yrand[i])*(y_sub-yrand[i]))
        
        # sky estimation
        sky_idx = np.logical_and(r>rap, r<rap + skywidth)
        sky = np.median(im_sub[sky_idx])
    
        # encircled flux
        phot_idx = np.where(r<=rap)
        total_pix = len(im_sub[phot_idx])
        flux = np.sum(im_sub[phot_idx]) - total_pix * sky 
        flux_arr.append(flux)

    # histogram 
    if hz1 == None:
        hz1 = min(flux_arr)
    if hz2 == None:
        hz2 = max(flux_arr)
    n, bins = np.histogram(np.array(flux_arr), bins=nhist, range=(hz1,hz2), normed=False, weights=None)
    
    # gaussian fit 
    bc = (bins[:-1] + bins[1:])/2
    p0 = [np.nanmax(n), 0.0, (hz2-hz1)/2.0]
    coeff, var_matrix = curve_fit(gauss, bc, n, p0=p0)
    n_fit = gauss(bc, *coeff)
    sigma = math.fabs(coeff[2])
    
    # plot 
    plt.ion()
    plt.figure()
    plt.step(bins[:-1], n, where='post')
    plt.plot(bc, n_fit)
    plt.xlabel('Aperture flux [ADU/sec]')
    plt.ylabel('Number')
    plt.title('Radius=%.1f, Nap=%d, Nbin=%d, Sigma=%.2f' % (rap, apnum, nhist, sigma))
    plt.grid()
    plt.show()
    _ = raw_input("Hit enter to finish the script:")        


if __name__=="__main__":

    usage = "usage: %prog input_image aperture_radius(pix) [options]"

    parser = OptionParser(usage)    
    parser.add_option("--xmin", dest="xmin", type="int", default=None,
                      help="Minimum X coordinate")
    parser.add_option("--xmax", dest="xmax", type="int", default=None,
                      help="Maximum X coordinate")
    parser.add_option("--ymin", dest="ymin", type="int", default=None,
                      help="Minimum Y coordinate")
    parser.add_option("--ymax", dest="ymax", type="int", default=None,
                      help="Maximum Y coordinate")
    parser.add_option("--skywidth", dest="skywidth", type="int", default=5,
                      help="Width of sky annulus in pix (default=5)")
    parser.add_option("--maxap", dest="maxap", type="int", default=500,
                      help="Maximum number of apertures (default=500)")
    parser.add_option("--nhist", dest="nhist", type="int", default=30,
                      help="Number of bins in the histogram (default=30)")
    parser.add_option("--hz1", dest="hz1", type="float", default=None,
                      help="Minimum bin value in the histogram")
    parser.add_option("--hz2", dest="hz2", type="float", default=None,
                      help="Maximum bin value in the histogram")
    
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()


    calcnoise(args[0], float(args[1]), xmin=options.xmin, xmax=options.xmax, ymin=options.ymin, ymax=options.ymax, skywidth=options.skywidth, maxap=options.maxap, nhist=options.nhist, hz1=options.hz1, hz2=options.hz2)

 
    
    
