#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
from optparse import OptionParser

def zcombine(inlst, outimg,combine='average',reject='none',masktype='none', review=False):

    imlist = []
    midpt = []

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
    ave = 0
    iraf.unlearn('imstat')
    for line in f:
        if not line.startswith('#'):
            fname = line[:-1]
            if os.access(fname, os.R_OK):
                imlist.append(fname)
                res = iraf.imstat(fname, format='no', fields='midpt', nclip=50, Stdout=1)
                if len(res) == 1:
                    midpt.append(float(res[0]))
                    ave += float(res[0])
                else:
                    print >> sys.stderr, 'failed to calculate midpt from the frame %s' % fname
                    f.close()
                    return 1
            else:
                print >> sys.stderr, 'cannot open image %s' % fname
                f.close()
                return 1
            nframe += 1

    f.close()

    if nframe > 0:
        ave = ave / nframe
    else:
        print >> sys.stderr, 'input list does not contain any frame' 
        return 1 

    # calculate zero shift value and record it into header
    iraf.unlearn('hedit')
    for i in range(nframe):
        midpt[i] = ave - midpt[i]
        iraf.hedit(imlist[i],'ZCOMBINE',midpt[i],add='yes',verify='no', show='no')
        
    # set parameters for imcombine 
    iraf.unlearn('imcombine')
    iraf.imcombine.combine = combine
    iraf.imcombine.reject = reject
    iraf.imcombine.zero = '!ZCOMBINE'
    iraf.imcombine.masktype = masktype

    # review imcombine parameters 
    if review:
        iraf.epar('imcombine')
    
    # execute imcombine
    iraf.imcombine('@'+inlst,outimg)

    # delete temporary added fits header 
    iraf.unlearn('hedit')
    for i in range(nframe):
        iraf.hedit(imlist[i],'ZCOMBINE','',delete='yes',verify='no', show='no')

    return 0

if __name__=="__main__":

    usage = "usage: %prog input_list output_image [options]"

    parser = OptionParser(usage)
    parser.add_option("--combine", dest="combine", type="choice", default="average",
                      choices=["average","median","sum"],
                      help="type of combine operation (average|median|sum)")
    parser.add_option("--reject", dest="reject", type="choice", default="none",
                      choices=["none", "minmax", "ccdclip", "crreject", "sigclip", "avsigclip", "pclipor"],
                      help="type of rejection (none|minmax|ccdclip|crreject|sigclip|avsigclip|pclipor)")
    parser.add_option("--masktype", dest="masktype", type="string", default="none",
                      help="mask type for imcombine (default=none)")
    parser.add_option("--review", dest="review", action="store_true", default=False,
                      help="review imcombine parameters (default=False)")
    
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()

    # execute 
    zcombine(args[0],args[1],combine=options.combine, reject=options.reject, masktype=options.masktype, review=options.review)
    
    

        
