#!/usr/bin/env python

import os,sys
import pyfits
from pyraf import iraf
from optparse import OptionParser

def pyhselect(img, hdrlist, lst):

    #<--- count number of lines in the list and backup the orignal file --->

    # open image list
    if lst:
        try :
            fimg=open(img)
        except IOError:
            print >> sys.stderr, 'cannot open %s' % img
            return 1
    
    # open header list
    try :
        fhdr=open(hdrlist)
    except IOError:
        print >> sys.stderr, 'cannot open %s' % hdrlist
        return 1

    # read header list
    hdr = '$I'
    for line in fhdr:
        if not line.startswith('#'):
            param = line.split()
            hdr = hdr + ',' + param[0]
    fhdr.close()

    if lst:
        for line in fimg:
            if not line.startswith('#'):
                im_file = line[:-1]
                iraf.hselect(im_file,hdr,expr='yes')
        fimg.close()
    else:
        #print hdr
        iraf.hselect(img,hdr,expr='yes')
        
    return 0


if __name__=="__main__":

    usage = "usage: %prog image header_list [options]"
    
    parser = OptionParser(usage)
    parser.add_option("--list", dest="lst", action="store_true", default=False,
                      help="input list instead of image (default=False)")

    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()
            
    pyhselect(args[0],args[1],lst=options.lst)
