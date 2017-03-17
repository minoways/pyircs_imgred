#!/usr/bin/env python
# IRCS imaging reduction pipeline (ver0.6) 
# Written by Y. Minowa 
# (( Revision History )) 
#  2013.10.10 -- Initial version (ver0.1)
#  2014.01.10 -- Revised after correcting bugs found in the quick reduction at the summit (ver0.2)
#  2014.02.12 -- Added erase option (ver0.3)
#  2015.10.01 -- Added fringe subtaction (ver.0.6)
import os,sys
from shutil import move
import pyfits
from pyraf import iraf
import tempfile 
from misc import *
from flatfield import * 
from distcor import *
from fringe_sub import *
from skysub_mod import *
from imshiftcomb import *
from optparse import OptionParser
from erase_file import * 
from copy_gmp import * 
from sigmap import *

def imgred_all_nb(inlist, output, start=0, end=15,
               flat='none', bpm='none', iternum=50, nsigma=3,
               dark='none',
               dist='none',
               fr_x0=0.0, fr_y0=0.0, fr_step=5,
               fitgeom='shift',
                  combine='average', reject='none', expmap='none', sigmap='none', whtmap='none', erase=False):

    ret = 0

    inpr = 'r'

    # erase all intermediate files 
    if erase:
        # step 1
        inpr = 'f2' + inpr
        if start < 2 and end >= 1:
            erase_file(inlist, inpref=inpr)

        # step 2
        inpr = 'fs' + inpr
        if start < 3 and end >= 2:
            erase_file(inlist, inpref=inpr)

        # step 3
        inpr = 's' + inpr
        if start < 4 and end >= 3:
            erase_file(inlist, inpref=inpr)

        # step 4
        if dist != 'none':
            inpr = 'd' + inpr
        if start < 5 and end >= 4:
            if dist != 'none':
                erase_file(inlist, inpref=inpr)

        # step 5
        if start < 6 and end >= 5:
            erase_file(inlist, inpref=inpr, ext='gmp')
            erase_file(inlist, inpref=inpr, ext='gmp2')
            erase_file(inlist, inpref=inpr, ext='dbs')
            if os.access(output, os.R_OK):
                os.remove(output)
            if expmap != 'none':
                if os.access(expmap, os.R_OK):
                    os.remove(expmap)
            if sigmap != 'none':
                if os.access(sigmap, os.R_OK):
                    os.remove(sigmap)

        return 0

    # 1. flat field
    if start < 2 and end >= 1:
        print '\n### Step 1: executing flat field ###\n'
        ret = flatfield(inlist, flat, inpref=inpr, outpref='f2', dark=dark, bpm=bpm)
        if ret != 0:
            print >> sys.stderr, 'Error in step 1'
            print >> sys.stderr, 'failed to execute flat field process' 
            return 1

    inpr = 'f2' + inpr 

    # 2. fringe subtraction
    if start < 3 and end >= 2:
        print '\n### Step 2: fringe subtraction ###\n'
        ret = fringe_sub(inlist, inpref=inpr, outpref='fs', mskpref=':OBJMASK', second=True, iternum=iternum, nsigma=nsigma, fr_x0=fr_x0, fr_y0=fr_y0, fr_step=fr_step)
        if ret != 0:
            print >> sys.stderr, 'Error in step 2'
            print >> sys.stderr, 'failed to subtract fringe from the flat corrected frames'
            return 1
            
    inpr = 'fs' + inpr

    # 3. subtract sky
    if start < 4 and end >= 3:
        print '\n### Step 3: subtracting sky ###\n'
        ret = skysub_mod(inlist, inpref=inpr, outpref='s', mskpref=':OBJMASK', iternum=iternum, nsigma=nsigma, step=0.05, second=True)
        if ret != 0:
            print >> sys.stderr, 'Error in step 3'
            print >> sys.stderr, 'failed to subtract sky from the fringe subtracted frames'
            return 1

    inpr = 's' + inpr
        
    # 4. distortion correction
    if start < 5 and end >= 4:
        if dist == 'none':
            print '\n### Skip Step 4: correcting distortion ###\n'
        else:
            print '\n### Step 4: correcting distortion for updated flat corrected frames###\n'
            ret = distcor(inlist, dist, inpref=inpr, outpref='d')
            if ret != 0:
                print >> sys.stderr, 'Error in step 4'
                print >> sys.stderr, 'failed to correct distortion'
                return 1

    if dist != 'none':
        inpr = 'd' + inpr

    
    # 5. shift and combine
    if start < 6 and end >= 5:
        ret = copy_gmp(inlist, inpref='sdfr', outpref=inpr)
        if ret != 0:
            print >> sys.stderr, 'Error in step 5'
            print >> sys.stderr, 'failed to copy gmp files'
            return 1                

        if os.access(output, os.R_OK):
            os.remove(output)
        print '\n### Step 5: combining frames (final step) ###\n'
        ret = imshiftcomb(inlist, output, fitgeom=fitgeom, inpref=inpr, objmask=':OBJMASK', combine=combine, reject=reject, second=True, first_pref='sdfr', second_pref=inpr)
        if ret != 0:
            print >> sys.stderr, 'Error in step 5'
            print >> sys.stderr, 'failed to shift and combine frames'
            return 1
        if sigmap != 'none':
            ret = sigmap(inlist, sigmap, expmap=expmap, whtmap=whtmap, inpref=inpr, ffpref='f2r', objmask=':OBJMASK', reject=reject)
            if ret != 0:
                print >> sys.stderr, 'Error in step 5'
                print >> sys.

    return 0
   

if __name__=="__main__":

    usage = "usage: %prog input_list output_image [options]"

    parser = OptionParser(usage)
    parser.add_option("--start", dest="start", type="int", default=0,
                      help="First step of the reduction sequence (default=0)")
    parser.add_option("--end", dest="end", type="int", default=15,
                      help="Last step of the reduction sequence (default=15)")
    parser.add_option("--flat", dest="flat", type="string", default="none",
                      help="Name of flat frame (default=none)")
    parser.add_option("--dark", dest="dark", type="string", default="none",
                      help="Name of dark frame (default=none)")
    parser.add_option("--dist", dest="dist", type="string", default="none",
                      help="Name of distortion database (default=none)")
    parser.add_option("--fitgeom", dest="fitgeom", type="choice", default="shift",
                      choices=["shift", "rotate"],
                      help="Fitting geometry for image shift (shift|rotate, default=shift)")
    parser.add_option("--fr_x0", dest="fr_x0", type="float", default=0.0,
                      help="Fringe center X coordinate (default=0.0)")
    parser.add_option("--fr_y0", dest="fr_y0", type="float", default=0.0,
                      help="Fringe center Y coordinate (default=0.0)")
    parser.add_option("--fr_step", dest="fr_step", type="int", default=5,
                      help="Step size for fringe fitting (default=5)")
    parser.add_option("--bpm", dest="bpm", type="string", default="none",
                      help="bad pixel mask frame (default=none)")
    parser.add_option("--iternum", dest="iternum", type="int", default=50,
                      help="Number of clipping iterations (default=50)")
    parser.add_option("--nsigma", dest="nsigma", type="float", default=3.0,
                      help="N-sigma reject limit (default=3.0)")
    parser.add_option("--combine", dest="combine", type="choice", default="average",
                      choices=["average","median","sum"],
                      help="type of combine operation (average|median|sum)")
    parser.add_option("--reject", dest="reject", type="choice", default="none",
                      choices=["none", "minmax", "ccdclip", "crreject", "sigclip", "avsigclip", "pclipor"],
                      help="type of rejection (none|minmax|ccdclip|crreject|sigclip|avsigclip|pclipor)")
    parser.add_option("--expmap", dest="expmap", type="string", default="none",
                      help="exposure map file name (default=none)")
    parser.add_option("--sigmap", dest="sigmap", type="string", default="none",
                      help="sigma map file name (default=none)")
    parser.add_option("--whtmap", dest="whtmap", type="string", default="none",
                      help="weight map file name (default=none)")
    parser.add_option("--erase", dest="erase", action="store_true", default=False,
                      help="erase intermediate files (default=False)")
    parser.add_option("--param", dest="param", type="string", default="none",
                      help="parameter file (default=none)")

    # parse options
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()

    # read parameter file
    if options.param != "none":
        if not os.path.exists(options.param):
            print >> sys.stderr, "parameter file (%s) does not exist" % (options.param)
            sys.exit()
        
        fparam = open(options.param)
        nline = 1
        for line in fparam:
            if not line.startswith('#'):
                pparam = line[:-1].split()
            
                if len(pparam) >= 2:
                    option_param = "--"+pparam[0]
                    if ''.join(sys.argv).find(option_param) < 0:
                        add_item = "--%s=%s" % (pparam[0],pparam[1])
                        sys.argv.append(add_item)
                else:
                    print >> sys.stderr, "Line %d: \"%s\" is ignored" % (nline, line[:-1])
            nline += 1
        fparam.close()
        
    # parse options again
    (options, args) = parser.parse_args()
        
    imgred_all_nb(args[0], args[1], start=options.start, end=options.end,
               flat=options.flat, bpm=options.bpm,
               iternum=options.iternum, nsigma=options.nsigma,
               dark=options.dark, dist=options.dist, 
               fr_x0=options.fr_x0, fr_y0=options.fr_y0, fr_step=options.fr_step,
               fitgeom=options.fitgeom,
                  combine=options.combine, reject=options.reject, expmap=options.expmap, sigmap=options.sigmap, whtmap=options.whtmap, erase=options.erase)
    
