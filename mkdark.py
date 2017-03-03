#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf

def mkdark(inlst, outimg):
    
    imlist = []
    exptime = []
    coadds = []
    ndr = []
    slwcnt = []

    # check if output file exists 
    if os.access(outimg, os.R_OK):
        print >> sys.stderr, 'operation would overwrite existing image (%s)' % outimg
        return 1

    # open input list and check if it exists
    try :
        f = open(inlst)
    except IOError:
        print >> sys.stderr, 'cannot open %s' % inlst
        return 1
    
    # check if file exist and record exptime, coadds, ndr, slwcnt
    nframe = 0
    for line in f:
        if not line.startswith('#'):
            fname = line[:-1]
            if os.access(fname, os.R_OK):
                imlist.append(fname)
                im = pyfits.open(fname)
                im_exptime = float(im[0].header['EXP1TIME'])
                im_coadds = int(im[0].header['COADDS'])
                im_ndr = int(im[0].header['NDR'])
                im_slwcnt = int(im[0].header['I_SLWCNT'])
                exptime.append(im_exptime)
                coadds.append(im_coadds)
                ndr.append(im_ndr)
                slwcnt.append(im_slwcnt)
                im.close()
            else:
                print >> sys.stderr, 'cannot open image %s' % fname
                f.close()
            nframe += 1
    f.close()

    # check consistency of EXPTIME, COADDS, NDR and SLWCNT thoughout frames in the list
    for i in range(nframe):
        if exptime[i] != exptime[0]:
            print >> sys.stderr, 'exptime for %s is not consistent with the first frame (%s)' % (imlist[i],imlist[0])
            return 1
        if coadds[i] != coadds[0]:
            print >> sys.stderr, 'coadds for %s is not consistent with the first frame (%s)' % (imlist[i],imlist[0])
            return 1
        if ndr[i] != ndr[0]:
            print >> sys.stderr, 'ndr for %s is not consistent with the first frame (%s)' % (imlist[i],imlist[0])
            return 1
        if slwcnt[i] != slwcnt[0]:
            print >> sys.stderr, 'slwcnt for %s is not consistent with the first frame (%s)' % (imlist[i],imlist[0])
            return 1
    
    # combine frames 
    iraf.unlearn('imcombine')
    iraf.imcombine.combine = 'median'
    iraf.imcombine.reject = 'sigclip'
    iraf.imcombine('@'+inlst, outimg)
    
    # divided by coadds and ndr
    iraf.unlearn('imarith')
    iraf.imarith(outimg,"/",ndr[0]*coadds[0],outimg,verbose='yes')

    return 0

if __name__=="__main__":

    argvs=sys.argv
    argc=len(argvs)

    if argc < 3:
        print 'Usage: %s [file_lst] [output image]' %(os.path.basename(sys.argv[0]))
        sys.exit()

    mkdark(sys.argv[1],sys.argv[2])
