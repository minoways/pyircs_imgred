#!/usr/bin/env python
import os,sys
import glob
import pyfits
import math
import numpy as np


def isfloat (instr):
    chkstr = str(instr)

    if chkstr[:1] == "-" or chkstr[:1] == "+":
        chkstr = chkstr[1:]
    numsplit = chkstr.split('.')
    if len(numsplit) == 1:
        if numsplit[0].isdigit():
            return True
        else:
            return False
    elif len(numsplit) == 2:
        if numsplit[0] == '':
            numsplit[0] = '0'
        if numsplit[1] == '':
            numsplit[1] = '0'
        if numsplit[0].isdigit() and numsplit[1].isdigit():
            return True
        else:
            return False
    else:
        return False

def remove_temp_all(prefix):
    files = glob.glob(prefix+'*')
    if len(files) > 0:
        for file in files:
            os.remove(file)


def check_input(inlist, inpref):

    # open input list and check if it exists
    try :
        fin = open(inlist)
    except IOError:
        print >> sys.stderr, 'cannot open %s' % inlist
        return 1

    inimg_arr = []
    for line in fin:
        if not line.startswith('#'):
            fname = line[:-1]
            if inpref == '':
                inimg = fname
            else:
                inimg = inpref + os.path.basename(fname)
            if not os.access(inimg, os.R_OK):
                print >> sys.stderr, 'cannot open image %s' % inimg
                fin.close()
                return 1 
            else:
                inimg_arr.append(inimg)
    fin.close()

    if len(inimg_arr) == 0:
        print >> sys.stderr, 'input list (%s) does not include any file' % inlist
        return 1

    return inimg_arr

def check_input_file(inlist, inpref, inext):

    # open input list and check if it exists
    try :
        fin = open(inlist)
    except IOError:
        print >> sys.stderr, 'cannot open %s' % inlist
        return 1

    infile_arr = []
    for line in fin:
        if not line.startswith('#'):
            fname, ext = os.path.splitext(os.path.basename(line[:-1]))
            infile = inpref + fname + inext
            if not os.access(infile, os.R_OK):
                print >> sys.stderr, 'cannot open file %s' % infile
                fin.close()
                return 1 
            else:
                infile_arr.append(infile)
    fin.close()

    if len(infile_arr) == 0:
        print >> sys.stderr, 'input list (%s) does not include any file' % inlist
        return 1

    return infile_arr

def check_inpref(inpref, inimg_arr):
    # check output prefix
    img_arr = []
    if inpref[0] == '@':
        try :
            f = open(inpref[1:])
        except IOError:
            print >> sys.stderr, 'cannot open input list (%s)' % inpref[1:]
            f.close()
            return 1

        for line in f:
            if not line.startswith('#'):
                img = line[:-1]
                if not os.access(img,os.R_OK):
                    print >> sys.stderr, 'cannot open image (%s)' % img
                    f.close()
                    return 1
                else :    
                    img_arr.append(img)
        f.close()

        if len(inimg_arr) != len(img_arr):
            print >> sys.stderr, 'number of images do not match with input image list'
            return 1
    elif inpref[0] == ':':
        for i in range(len(inimg_arr)):
            im = pyfits.open(inimg_arr[i])
            try:
                img = im[0].header[inpref[1:]]
            except KeyError:
                print sys.stderr, 'cannot read header parameter (%s)' % inpref[1:]
                im.close()
                return 1
            im.close()
            imgfile = img.split('[')
            if not os.access(imgfile[0],os.R_OK):
                print >> sys.stderr, 'cannot open image (%s)' % img
                return 1
            else :    
                img_arr.append(imgfile[0])
    else :
        for i in range(len(inimg_arr)):
            img = inpref + os.path.basename(inimg_arr[i])
            if not os.access(img,os.R_OK):
                print >> sys.stderr, 'cannot open image (%s)' % img
                return 1
            else :    
                img_arr.append(img)
    
    return img_arr

def check_outpref(outpref, inimg_arr):
    # check output prefix
    outimg_arr = []
    if outpref[0] == '@':
        try :
            fout = open(outpref[1:])
        except IOError:
            print >> sys.stderr, 'cannot open output list (%s)' % outpref[1:]
            fout.close()
            return 1

        for line in fout:
            if not line.startswith('#'):
                outimg = line[:-1]
                if os.access(outimg,os.R_OK):
                    print >> sys.stderr, 'operation would overwrite exisiting image (%s)' % outimg
                    fout.close()
                    return 1
                else :    
                    outimg_arr.append(outimg)
        fout.close()

        if len(inimg_arr) != len(outimg_arr):
            print >> sys.stderr, 'number of input and output images do not match'
            return 1
    else :
        for i in range(len(inimg_arr)):
            outimg = outpref + os.path.basename(inimg_arr[i])
            if os.access(outimg,os.R_OK):
                print >> sys.stderr, 'operation would overwrite exisiting image (%s)' % outimg
                return 1
            else :    
                outimg_arr.append(outimg)
    
    return outimg_arr

def robust_sigma(in_y, zero=0):
   """
   Calculate a resistant estimate of the dispersion of
   a distribution. For an uncontaminated distribution,
   this is identical to the standard deviation.
 
   Use the median absolute deviation as the initial
   estimate, then weight points using Tukey Biweight.
   See, for example, Understanding Robust and
   Exploratory Data Analysis, by Hoaglin, Mosteller
   and Tukey, John Wiley and Sons, 1983.
 
   .. note:: ROBUST_SIGMA routine from IDL ASTROLIB.
 
   :History:
       * H Freudenreich, STX, 8/90
       * Replace MED call with MEDIAN(/EVEN), W. Landsman, December 2001
       * Converted to Python by P. L. Lim, 11/2009
 
   Examples
   --------
   >>> result = robust_sigma(in_y, zero=1)
 
   Parameters
   ----------
   in_y: array_like
       Vector of quantity for which the dispersion is
       to be calculated
 
   zero: int
       If set, the dispersion is calculated w.r.t. 0.0
       rather than the central value of the vector. If
       Y is a vector of residuals, this should be set.
 
   Returns
   -------
   out_val: float
       Dispersion value. If failed, returns -1.
 
   """
   # Flatten array
   y = in_y.reshape(in_y.size, )
 
   eps = 1.0E-20
   c1 = 0.6745
   c2 = 0.80
   c3 = 6.0
   c4 = 5.0
   c_err = -1.0
   min_points = 3
 
   if zero:
       y0 = 0.0
   else:
       y0 = np.median(y)
 
   dy    = y - y0
   del_y = abs( dy )
 
   # First, the median absolute deviation MAD about the median:
 
   mad = np.median( del_y ) / c1
 
   # If the MAD=0, try the MEAN absolute deviation:
   if mad < eps:
       mad = np.mean( del_y ) / c2
   if mad < eps:
       return 0.0
 
   # Now the biweighted value:
   u  = dy / (c3 * mad)
   uu = u*u
   q  = np.where(uu <= 1.0)
   count = len(q[0])
   if count < min_points:
       print 'ROBUST_SIGMA: This distribution is TOO WEIRD! Returning', c_err
       return c_err
 
   numerator = np.sum( (y[q]-y0)**2.0 * (1.0-uu[q])**4.0 )
   n    = y.size
   den1 = np.sum( (1.0-uu[q]) * (1.0-c4*uu[q]) )
   siggma = n * numerator / ( den1 * (den1 - 1.0) )
 
   if siggma > 0:
       out_val = np.sqrt( siggma )
   else:
       out_val = 0.0
 
   return out_val

def meanclip(indata, clipsig=3.0, maxiter=5, converge_num=0.001, median=0, verbose=0):
   """
   Computes an iteratively sigma-clipped mean on a
   data set. Clipping is done about median, but mean
   is returned.
 
   .. note:: MYMEANCLIP routine from ACS library.
 
   :History:
       * 21/10/1998 Written by RSH, RITSS
       * 20/01/1999 Added SUBS, fixed misplaced paren on float call, improved doc. RSH
       * 24/11/2009 Converted to Python. PLL.
 
   Examples
   --------
   >>> mean, sigma = meanclip(indata)
 
   Parameters
   ----------
   indata: array_like
       Input data.
 
   clipsig: float
       Number of sigma at which to clip.
 
   maxiter: int
       Ceiling on number of clipping iterations.
 
   converge_num: float
       If the proportion of rejected pixels is less than
       this fraction, the iterations stop.
 
   median: {0, 1}
       Calculate median instead of mean?

   verbose: {0, 1}
       Print messages to screen?
 
   Returns
   -------
   mean: float
       N-sigma clipped mean.
 
   sigma: float
       Standard deviation of remaining pixels.
 
   """
   # Flatten array
   skpix = indata.reshape( indata.size, )

   # mask nan
   skpix = np.ma.masked_array(skpix, np.isnan(skpix))

   ct = indata.size
   iter = 0; c1 = 1.0 ; c2 = 0.0
 
   while (c1 >= c2) and (iter < maxiter):
       lastct = ct
       medval = np.median(skpix)
       sig = np.std(skpix)
       #print '%f %f %f' % (clipsig, sig, medval)
       wsm = np.where( abs(skpix-medval) < clipsig*sig )
       ct = len(wsm[0])
       if ct > 0:
           skpix = skpix[wsm]
 
       c1 = abs(ct - lastct)
       c2 = converge_num * lastct
       iter += 1
   # End of while loop
 
   if median:
       mean  = np.median( skpix )
   else:
       mean  = np.mean( skpix )
   sigma = robust_sigma( skpix )
 
   if verbose:
       prf = 'MEANCLIP:'
       if median:
           print '%s %.1f-sigma clipped median' % (prf, clipsig)
           print '%s Median computed in %i iterations' % (prf, iter)
           print '%s Median = %.6f, sigma = %.6f' % (prf, mean, sigma)
       else:
           print '%s %.1f-sigma clipped mean' % (prf, clipsig)
           print '%s Mean computed in %i iterations' % (prf, iter)
           print '%s Mean = %.6f, sigma = %.6f' % (prf, mean, sigma)
   return mean, sigma

def convert_maskfits_int(input, output):
    # open fits file
    im = pyfits.open(input)
    # convert the data to nearest int 
    data = np.rint(im[0].data)
    # close image 
    im.close()
    
    # output fits 
    hdu = pyfits.PrimaryHDU(data)
    out = pyfits.HDUList([hdu])
    if input == output:
        os.remove(input)
    out.writeto(output)


