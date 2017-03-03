#!/usr/bin/env python 
import os,sys
from pyraf import iraf
import numpy as np
import tempfile
from optparse import OptionParser
from ds9 import ds9 
from misc import *
import math 

def galfit_psf(image, output):
    rap = 20
    skycal = "const"
    annu = 20
    dannu = 5
    mask = "none"
    constraint = "none"
    sigma = "none"
    sr_init = 0.3

    # check input image
    if not os.access(image, os.R_OK):
        print >> sys.stderr, 'cannot read input image (%s)' % image
        return 1

    # check output
    if os.access(output, os.R_OK):
        print >> sys.stderr, 'operation would overwrite existing data (%s)' % output
        return 1
    else:
        fout = open(output,'w')
    
    # temporary file
    tmp = tempfile.NamedTemporaryFile(suffix='', prefix='', dir='/tmp')
    tmp_prefix = tmp.name
    tmp.close()

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
    iraf.display(image, 1)

    # imexam
    print '\nSelect Object --> type a , Quit --> type q on the ds9 image'
    print 'Caution: do not type any key except a or q'
    ret = iraf.imexam(image, 1, Stdout=1)
    
    # display result
    print '### Imexamine ###'
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
            print >> sys.stderr, 'failed to pick up object in %s' % image
            return 1
    
    # read data into memory
    im = pyfits.open(image)
    im_data = im[0].data
    im_hdr = im[0].header
    im.close()

    # coordinates and radius
    x = np.ones(im_data.shape[0]*im_data.shape[1]).reshape(im_data.shape[0],im_data.shape[1])* np.arange(im_data.shape[1]) + 1
    y = (np.ones(im_data.shape[1]*im_data.shape[0]).reshape(im_data.shape[1],im_data.shape[0]) * np.arange(im_data.shape[0]) + 1).T
    r = np.sqrt((x-xobj)*(x-xobj) + (y-yobj)*(y-yobj))

    # sky estimation
    if skycal == "constant":
        sky = skyconst
    elif skycal == "median":
        sky_idx = np.logical_and(r>=annu, r<annu+dannu)
        sky = np.median(im_data[sky_idx])
    else:
        sky = 0.0

    # aperture photometry
    phot_idx = np.where(r<=rap)
    total_pix = len(im_data[phot_idx])
    flux_total = np.sum(im_data[phot_idx]) - total_pix * sky 
    print '### Aperture Photometry ###'
    print "total_flux = %f [ADU], sky_level = %f [ADU/pux]" % (flux_total, sky) 
  
    flux_core_init = sr_init * flux_total
    flux_halo_init = (1.0 - sr_init) * flux_total
    mag_core_init = 25.0 - 2.5*math.log10(flux_core_init)
    mag_halo_init = 25.0 - 2.5*math.log10(flux_halo_init)

    # check sampling
    mag_fac = 1
    fwhm_org = fwhm
    while fwhm < 3.0:
        mag_fac = mag_fac * 2
        fwhm = fwhm_org * mag_fac

    xobj = xobj * mag_fac
    yobj = yobj * mag_fac

    # magnify image
    if mag_fac > 1:
        tmp_resample = '%s_rs.fits' % (tmp_prefix)
        iraf.magnify(image, tmp_resample, mag_fac, mag_fac, interpo="drizzle", boundar="constant", constant=0, fluxcon="yes")
    else:
        tmp_resample = image

    # create sigma image
    tmp_sigma = '%s_sigma.fits' % (tmp_prefix)
    if sigma == 'none':
        gain = im_hdr['Gain']
        iraf.imarith(tmp_resample, "*", gain, tmp_sigma)
        iraf.stsdas.toolbox.imgtools
        iraf.imcalc(tmp_sigma, tmp_sigma, "sqrt(im1)")
        iraf.imarith(tmp_sigma, "/", gain, tmp_sigma)
    else:
        iraf.imarith(sigma,"*",sigma, tmp_sigma)
        iraf.magnify(tmp_sigma, tmp_sigma, mag_fac, mag_fac, interpo="drizzle", boundar="constant", constant=0, fluxcon="yes")
        iraf.imcalc(tmp_sigma, tmp_sigma, "sqrt(im1)")

    # fitting control parameters 
    tmp_galfit = '%s_galfit.txt' % (tmp_prefix)
    tmp_imgblock = '%s_imgblock.fits' % (tmp_prefix)
    f_galfit = open(tmp_galfit,'w')
    f_galfit.write('A) %s\n' % tmp_resample)
    f_galfit.write('B) %s\n' % tmp_imgblock)
    f_galfit.write('C) %s\n' % tmp_sigma)
    f_galfit.write('D) none\n')
    f_galfit.write('E) 1\n')
    if mask != 'none':
        f_galfit.write('F) %s\n' % mask)
    else:
        f_galfit.write('F) none\n')
    if constraint != 'none':
        f_galfit.write('G) %s\n' % constraint)
    else:
        f_galfit.write('G) none\n')
    f_galfit.write('H) %d %d %d %d\n' % (xmin,xmax,ymin,ymax))
    f_galfit.write('I) 100 100\n' % tmp_resample)
    f_galfit.write('J) 25.0\n' % tmp_resample)
    f_galfit.write('K) %f %f\n' % (pixscale, pixscale))
    f_galfit.write('O) regular\n' % tmp_resample)
    f_galfit.write('P) 0\n' % tmp_resample)
    f_galfit.write('\n' % tmp_resample)
    
    # core fitting initial parameters 
    if core_prof == 'moffat':
        f_galfit.write('0) moffat\n')
    else:
        f_galfit.write('0) gaussian\n')
    f_galfit.write('2) %f %f 1 1\n' % (xobj, yobj))
    f_galfit.write('3) %f 1\n' % mag_core_init)
    f_galfit.write('4) %f 1\n' % fwhm)
    if core_prof == 'moffat':
        f_galfit.write('5) %f 1\n' % beta_init)
    else:
        f_galfit.write('5) 0.0000 0\n')
    f_galfit.write('6) 0.0000 0\n')
    f_galfit.write('7) 0.0000 0\n')
    f_galfit.write('8) 0.0000 0\n')
    f_galfit.write('9) 0.9 1\n')
    f_galfit.write('10) 0.0 1\n')
    f_galfit.write('Z) 0\n')

    # halo fitting initial parameters 
    f_galfit.write('0) moffat\n')
    f_galfit.write('2) %f %f 1 1\n' % (xobj, yobj))
    f_galfit.write('3) %f 1\n' % mag_halo_init)
    f_galfit.write('4) %f 1\n' % fwhm)
    f_galfit.write('5) %f 1\n' % beta_init)
    f_galfit.write('6) 0.0000 0\n')
    f_galfit.write('7) 0.0000 0\n')
    f_galfit.write('8) 0.0000 0\n')
    f_galfit.write('9) 0.9 1\n')
    f_galfit.write('10) 0.0 1\n')
    f_galfit.write('Z) 0\n')
    
    # sky fitting parameters 
    if skycal == "fit":
        f_galfit.write('0) sky\n')
        f_galfit.write('1) 0.000      1\n')
        f_galfit.write('2) 0.000      0\n')
        f_galfit.write('3) 0.000      0\n')
        f_galfit.write('Z) 0\n')

    f_galfit.close()
    
    # execute galfit
    cmd = 'galfit %s' % tmp_galfit
    subprocess.call(cmd, shell=True)

if __name__=="__main__":

    usage = "usage: %prog inimage output [options]"

    parser = OptionParser(usage)
    #parser.add_option("--inpref", dest="inpref", type="string", default="",
    #                  help="prefix for input frame (default=none)")
    #parser.add_option("--outpref", dest="outpref", type="string", default="s",
    #                  help="prefix for output frame or output list with @ mark (default=s)")
    #parser.add_option("--mskpref", dest="mskpref", type="string", default="none",
    #                  help="prefix for mask image, image list with \'@\', or header keyword with \':\'")
    #parser.add_option("--ndthpos", dest="ndthpos", type="int", default=9,
    #                  help="number of dithering position (default=9)")
    #parser.add_option("--iternum", dest="iternum", type="int", default=50,
    #                  help="Number of clipping iterations (default=50)")
    #parser.add_option("--nsigma", dest="nsigma", type="float", default=3.0,
    #                  help="N-sigma reject limit (default=3.0)")
    #parser.add_option("--second", dest="second", action="store_true", default=False,
    #                  help="Is this second cycle subtraction? (default=False)")
    #parser.add_option("--step", dest="step", type="float", default=0.05,
    #                  help="step for finding the best scale value (default=0.05)")
    #parser.add_option("--fringe", dest="fringe", action="store_true", default=False,
    #                  help="Remove fringe pattern? (default=False)")
    #parser.add_option("--fr_x0", dest="fr_x0", type="float", default=0.0,
    #                  help="Fringe center X coordinate (default=0.0)")
    #parser.add_option("--fr_y0", dest="fr_y0", type="float", default=0.0,
    #                  help="Fringe center Y coordinate (default=0.0)")
    #parser.add_option("--fr_step", dest="fr_step", type="int", default=5,
    #                  help="Step size for fringe fitting (default=5)")

    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()
                
    galfit_psf(args[0], args[1])
  
