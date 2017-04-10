#!/usr/bin/env python
# IRCS imaging reduction pipeline (ver0.6) 
# Written by Y. Minowa 
# (( Revision History )) 
#  2013.10.10 -- Initial version (ver0.1)
#  2014.01.10 -- Revised after correcting bugs found in the quick reduction at the summit (ver0.2)
#  2014.02.12 -- Added erase option (ver0.3)
#  2015.10.01 -- Added fringe subtaction (ver.0.6)
#  2017.03.08 -- Added parameter file option
import os,sys
import os.path
from shutil import move
import pyfits
from pyraf import iraf
import tempfile 
from misc import *
from mkskyflat import *
from mkskyflat_quick import *
from flatfield import * 
from mkobjmask import * 
from distcor import *
from skysub import *
from skysub_mod import *
from calcoffset import *
from imshiftcomb import *
from mkskyflat_with_mask import *
from registmask import *
from mkcombmask import *
from copy_frame import *
from optparse import OptionParser
from comment_out_skipped import *
from erase_file import * 
from copy_gmp import * 
from sigmap import *

def imgred_all(inlist, output, quick=False, start=0, end=15,
               flat='none', skyflat=True, bpm='none', iternum=50, nsigma=3, minpix=250, hsigma=3, lsigma=10, conv='block 3 3',
               dark='none',
               dist='none',
               fringe=False, fr_x0=0.0, fr_y0=0.0, fr_step=5,
               trace=True, fitgeom='shift',
               combine='average', reject='none', expmap='none', sigmap='none', whtmap='none', erase=False):
    
    ret = 0

    inpr = ''

    # erase all intermediate files 
    if erase:
        # step 0
        inpr = 'r' + inpr
        if start < 1 and end >= 0:
            erase_file(inlist, inpref=inpr)

        # step 1
        if start < 2 and end >= 1:
            erase_file(inlist, inpref='m'+inpr)
            if skyflat:
                if os.access(flat, os.R_OK):
                    iraf.delete(flat, verify='no')

        # step 2
        inpr = 'f' + inpr
        if start < 3 and end >= 2:
            erase_file(inlist, inpref=inpr)

        # step 3
        if dist != 'none':
            inpr = 'd' + inpr
        if start < 4 and end >= 3:
            if dist != 'none':
                erase_file(inlist, inpref=inpr)

        # step 4
        if start < 5 and end >= 4:
            erase_file(inlist, inpref='m'+inpr)

        # step 5
        inpr = 's' + inpr
        if start < 6 and end >= 5:
            erase_file(inlist, inpref='sdfr')

        # step 6
        if start < 7 and end >= 6:
            erase_file(inlist, inpref='sdfr', ext='gmp')
            
        # step 7
        if start < 8 and end >= 7:
            if os.access(inlist+".skiporg", os.R_OK):
                shutil.move(inlist+".skiporg",inlist)

        # step 8
        if start < 9 and end >= 8:
            erase_file(inlist, inpref='sdfr', ext='gmp2')
            erase_file(inlist, inpref='sdfr', ext='dbs')
            if os.access(output, os.R_OK):
                iraf.delete(output, verify='no')
            if sigmap != 'none' and os.access(sigmap, os.R_OK):
                iraf.delete(sigmap, verify='no')

        # step 9
        inpr = 'r'
        if start < 10 and end >= 9:
            if dist == 'none':
                erase_file(inlist, inpref='mcsf'+inpr)
            else:
                erase_file(inlist, inpref='mcsdf'+inpr)
                
        # step 10 (no intermediate file)
        if start < 11 and end >= 10:
            pass

        # step 11
        if start < 12 and end >= 11:
            if skyflat:
                if os.access(flat, os.R_OK):
                    iraf.delete(flat, verify='no')

        # step 12
        inpr = 'f2' + inpr
        if start < 13 and end >= 12:
            erase_file(inlist, inpref=inpr)

        # step 13
        if dist != 'none':
            inpr = 'd' + inpr
        if start < 14 and end >= 13:
            if dist != 'none':
                erase_file(inlist, inpref=inpr)
                
        # step 14 
        inpr = 's' + inpr
        if start < 15 and end >= 14:
            erase_file(inlist, inpref=inpr)

        # step 15
        if start < 16 and end >= 15:
            erase_file(inlist, inpref='sdf2r', ext='gmp')
            erase_file(inlist, inpref='sdf2r', ext='gmp2')
            erase_file(inlist, inpref='sdf2r', ext='dbs')
            if os.access(output, os.R_OK):
                os.remove(output)
            if expmap != 'none':
                if os.access(expmap, os.R_OK):
                    os.remove(expmap)
            if whtmap != 'none':
                if os.access(whtmap, os.R_OK):
                    os.remove(whtmap)
            if sigmap != 'none':
                if os.access(sigmap, os.R_OK):
                    os.remove(sigmap)

        return 0

    # 0. copy raw frames to the current directory
    if start < 1 and end >= 0:
        print '\n### Step 0: copy input frames to the current directory ###\n'
        ret = copy_frame(inlist, inpref=inpr, outpref='r')
        if ret != 0:
            print >> sys.stderr, 'Error in step 0'
            print >> sys.stderr, 'failed to copy input frames to the current directory'
            return 1
    inpr = 'r' + inpr

    # 1. making sky flat
    if start < 2 and end >= 1:
        print '\n### Step 1: making self sky flat frame ###\n'
        if skyflat:
            if quick:
                ret = mkskyflat_quick(inlist, flat, inpref=inpr, bpm=bpm, iternum=iternum, nsigma=nsigma)
                if ret != 0:
                    print >> sys.stderr, 'Error in step 1'
                    print >> sys.stderr, 'failed to make self sky flat frame'
                    return 1
            else:
                ret = mkskyflat(inlist, flat, inpref=inpr, mskpref='m', bpm=bpm, iternum=iternum, nsigma=nsigma, minpix=minpix, hsigma=hsigma, lsigma=lsigma, conv=conv)
                if ret != 0:
                    print >> sys.stderr, 'Error in step 1'
                    print >> sys.stderr, 'failed to make self sky flat frame'
                    return 1
        else:
            print 'Skip step 1 and use %s as a flat frame' % (flat)
   
    # 2. flat field
    if start < 3 and end >= 2:
        print '\n### Step 2: executing flat field ###\n'
        ret = flatfield(inlist, flat, inpref=inpr, outpref='f', dark=dark, bpm=bpm)
        if ret != 0:
            print >> sys.stderr, 'Error in step 2'
            print >> sys.stderr, 'failed to execute flat field process' 
            return 1

    inpr = 'f' + inpr 

    # 3. distortion correction
    if start < 4 and end >= 3:
        if dist == 'none':
            print '\n### Skip Step 3: correcting distortion ###\n'
        else:
            print '\n### Step 3: correcting distortion ###\n'
            ret = distcor(inlist, dist, inpref=inpr, outpref='d')
            if ret != 0:
                print >> sys.stderr, 'Error in step 3'
                print >> sys.stderr, 'failed to correct distortion'
                return 1

    if dist != 'none':
        inpr = 'd' + inpr

    # 4. make object mask 
    if start < 5 and end >= 4:
        print '\n### Step 4: making object mask ###\n'
        ret = mkobjmask(inlist, inpref=inpr, outpref='m', bpm=bpm, minpix=minpix, hsigma=hsigma, lsigma=lsigma, conv=conv)
        if ret != 0:
            print >> sys.stderr, 'Error in step 4'
            print >> sys.stderr, 'failed to make object mask for distortion corrected frames'
            return 1
   
    # 5. subtract sky
    if start < 6 and end >= 5:
        print '\n### Step 5: subtracting sky ###\n'
        if not quick:
            ret = skysub_mod(inlist, inpref=inpr, outpref='s', mskpref=':OBJMASK', iternum=iternum, nsigma=nsigma, step=0.05, second=False)
        else:
            ret = skysub(inlist, inpref=inpr, outpref='s', mskpref=':OBJMASK', iternum=iternum, nsigma=nsigma)
        #ret = skysub(inlist, inpref=inpr, outpref='s', mskpref=':OBJMASK', iternum=iternum, nsigma=nsigma)
        if ret != 0:
            print >> sys.stderr, 'Error in step 5'
            print >> sys.stderr, 'failed to subtract sky from the distortion corrected frames'
            return 1

    inpr = 's' + inpr

    # 6. measure offset
    if start < 7 and end >= 6:
        print '\n### Step 6: measuring offset ###\n'
        #ret = calcoffset(inlist, offset=offset, inpref=inpr, xc=xc, yc=yc, trace=trace, review=False)
        ret = calcoffset(inlist, inpref=inpr, trace=trace, review=False)
        if ret != 0:
            print >> sys.stderr, 'Error in step 6'
            print >> sys.stderr, 'failed to measure offset'
            return 1

    # 7. comment out skipped image
    if start < 8 and end >= 7:
        print '\n### Step 7: comment out skipped images ###\n'
        ret = comment_out_skipped(inlist, sspref=inpr)
        if ret != 0:
            print >> sys.stderr, 'Error in step 7'
            print >> sys.stderr, 'failed to comment out skipped line'
            return 1

    # 8. shift and combine
    if start < 9 and end >= 8:
        print '\n### Step 8: combining frames ###\n'
        ret = imshiftcomb(inlist, output, fitgeom=fitgeom, inpref=inpr, objmask=':OBJMASK', combine=combine, reject=reject, second=False, sigmap=sigmap, ffpref='fr')
        if ret != 0:
            print >> sys.stderr, 'Error in step 8'
            print >> sys.stderr, 'failed to shift and combine frames'
            return 1
    
    if not quick:
        # 9. make object mask from the combined image
        if start < 10 and end >= 9:
            print '\n### Step 9: making object mask from the combined image ###\n'
            ret = mkcombmask(output, inlist, sspref=inpr, outpref='mc', minpix=minpix, hsigma=hsigma, lsigma=lsigma, conv=conv, bpm='none', sigmap=sigmap)
            if ret != 0:
                print >> sys.stderr, 'Error in step 9'
                print >> sys.stderr, 'failed to make object mask from the combined image'
                return 1
        
        inpr = 'r'
        if dist == 'none':
            mskpref = 'mcsf'
        else:
            mskpref = 'mcsdf'

        # 10. register object mask into the fits header of raw frames
        if start < 11 and end >= 10:
            print '\n### Step 10: registering object mask frames into the fits header of raw frames ###\n'
            ret = registmask(inlist, inpref=inpr, mskpref=mskpref, key='OBJMASK')
            if ret != 0:
                print >> sys.stderr, 'Error in step 10'
                print >> sys.stderr, 'failed to register  object mask from the combined image'
                return 1

        # 11. make self sky flat frame
        if start < 12 and end >= 11:
            if skyflat:
                print '\n### Step 11: making self sky flat frame with updated object mask frames ###\n'
                if os.access(flat, os.R_OK):
                    os.remove(flat)
                ret = mkskyflat_with_mask(inlist, flat, inpref=inpr, mskpref=':OBJMASK', bpm=bpm, iternum=iternum, nsigma=nsigma)
                if ret != 0:
                    print >> sys.stderr, 'Error in step 11'
                    print >> sys.stderr, 'failed to register  object mask from the combined image'
                    return 1
            else:
                print 'Skip step 11 use %s as a flat frame' % (flat)
            
        # 12. correct flat field with the new flat frame
        if start < 13 and end >= 12:
            print '\n### Step 12: executing flat field with updated flat frame ###\n'
            ret = flatfield(inlist, flat, inpref=inpr, outpref='f2', dark=dark, bpm=bpm)
            if ret != 0:
                print >> sys.stderr, 'Error in step 12'
                print >> sys.stderr, 'failed to execute flat field process' 
                return 1

        inpr = 'f2' + inpr

        # 13. distortion correction
        if start < 14 and end >= 13:
            if dist == 'none':
                print '\n### Skip Step 13: correcting distortion for updated flat corrected frames###\n'
            else:
                print '\n### Step 13: correcting distortion for updated flat corrected frames###\n'
                ret = distcor(inlist, dist, inpref=inpr, outpref='d')
                if ret != 0:
                    print >> sys.stderr, 'Error in step 13'
                    print >> sys.stderr, 'failed to correct distortion'
                    return 1

        if dist != 'none':
            inpr = 'd' + inpr

        # 14. subtract sky 
        if start < 15 and end >= 14:
            print '\n### Step 14: subtracting sky ###\n'
            ret = skysub_mod(inlist, inpref=inpr, outpref='s', mskpref=':OBJMASK', iternum=iternum, nsigma=nsigma, step=0.05, second=True)
            if ret != 0:
                print >> sys.stderr, 'Error in step 14'
                print >> sys.stderr, 'failed to subtract sky from the distortion corrected frames'
                return 1

        inpr = 's' + inpr

        # 15. shift and combine
        if start < 16 and end >= 15:
            
            ret = copy_gmp(inlist, inpref='sdfr', outpref=inpr)
            if ret != 0:
                print >> sys.stderr, 'Error in step 15'
                print >> sys.stderr, 'failed to copy gmp files'
                return 1                

            if os.access(output, os.R_OK):
                os.remove(output)
            print '\n### Step 15: combining frames (final step) ###\n'
            ret = imshiftcomb(inlist, output, fitgeom=fitgeom, inpref=inpr, objmask=':OBJMASK', combine=combine, reject=reject, second=True, ffpref='f2r', expmap=expmap, whtmap=whtmap, sigmap=sigmap)
            if ret != 0:
                print >> sys.stderr, 'Error in step 15'
                print >> sys.stderr, 'failed to shift and combine frames'
                return 1
                
    return 0


if __name__=="__main__":

    usage = "usage: %prog input_list output_image [options]"

    parser = OptionParser(usage)
    parser.add_option("--start", dest="start", type="int", default=0,
                      help="First step of the reduction sequence (default=0)")
    parser.add_option("--end", dest="end", type="int", default=15,
                      help="Last step of the reduction sequence (default=15)")
    parser.add_option("--quick", dest="quick", action="store_true", default=False,
                      help="execute quick reduction (default=False)")
    parser.add_option("--flat", dest="flat", type="string", default="none",
                      help="Name of flat frame (default=none)")
    parser.add_option("--skyflat", dest="skyflat", action="store_true", default=False,
                      help="Make self sky flat frame (default=False)")
    parser.add_option("--dark", dest="dark", type="string", default="none",
                      help="Name of dark frame (default=none)")
    parser.add_option("--dist", dest="dist", type="string", default="none",
                      help="Name of distortion database (default=none)")
    parser.add_option("--fitgeom", dest="fitgeom", type="choice", default="shift",
                      choices=["shift", "rotate"],
                      help="Fitting geometry for image shift (shift|rotate, default=shift)")
    #parser.add_option("--offset", dest="offset", type="string", default="none",
    #                  help="offset data file (default=none)")
    #parser.add_option("--xc", dest="xc", type="float", default=512.0,
    #                  help="origin of X coordinate (default=512.0)")
    #parser.add_option("--yc", dest="yc", type="float", default=512.0,
    #                  help="origin of Y coordinate (default=512.0)")
    #parser.add_option("--ndthpos", dest="ndthpos", type="int", default=9,
    #                  help="number of dithering position (default=9)")
    parser.add_option("--fringe", dest="fringe", action="store_true", default=False,
                      help="Remove fringe pattern? (default=False)")
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
    parser.add_option("--minpix", dest="minpix", type="int", default=250,
                      help="minimum number of pixels in detected objects (default=250)")
    parser.add_option("--hsigma", dest="hsigma", type="float", default=3.0,
                      help="sigma threshold above sky (default=3)")
    parser.add_option("--lsigma", dest="lsigma", type="float", default=10.0,
                      help="sigma threshold below sky (default=10)")
    parser.add_option("--conv", dest="conv", type="string", default="block 3 3",
                      help="convolution kernel for objmask (default=\"block 3 3\")")
    parser.add_option("--combine", dest="combine", type="choice", default="average",
                      choices=["average","median","sum"],
                      help="type of combine operation (average|median|sum)")
    parser.add_option("--reject", dest="reject", type="choice", default="sigclip",
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
    parser.add_option("--notrace", dest="notrace", action="store_false", default=True,
                      help="automatic trace of object position (default=True)")
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

    # execute pipeline    
    imgred_all(args[0], args[1], quick=options.quick, start=options.start, end=options.end,
               flat=options.flat, skyflat=options.skyflat, bpm=options.bpm,
               iternum=options.iternum, nsigma=options.nsigma, minpix=options.minpix, hsigma=options.hsigma, lsigma=options.lsigma, conv=options.conv,
               dark=options.dark, dist=options.dist, 
               fringe=options.fringe,fr_x0=options.fr_x0, fr_y0=options.fr_y0, fr_step=options.fr_step,
               trace=options.notrace, fitgeom=options.fitgeom,
               combine=options.combine, reject=options.reject, expmap=options.expmap, sigmap=options.sigmap, whtmap=options.whtmap, erase=options.erase)
    
