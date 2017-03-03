#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
import tempfile 
from misc import *
from optparse import OptionParser

def distcor(inlist, database, inpref='', outpref='d'):
    
    # open input list and check if it exists
    inimg_arr = check_input(inlist, inpref)
    if isinstance(inimg_arr, int):
        return 1

    # check output prefix
    outimg_arr = check_outpref(outpref, inimg_arr)
    if isinstance(inimg_arr,int):
        return 1

    # check database 
    try :
        fd = open(database)
    except IOError:
        print >> sys.stderr, 'cannot open database file (%s)' % database
        f.close()
        return 1
    
    # read geomap coordinate file
    read = 0
    for line in fd:
        param = line[:-1].split()
        if param[0] == 'begin':
            gmp = param[1]
            read += 1
    fd.close()

    if read == 0:
        print >> sys.stderr, 'cannot read coordinate file in geomap database (%s)' % database
        f.close
        return 1

    # correct distortion 
    for i in range(len(inimg_arr)):
        # correct distortion
        iraf.unlearn('geotran')
        iraf.geotran(inimg_arr[i], outimg_arr[i], database, gmp)

    return 0

if __name__=="__main__":

    usage = "usage: %prog input_list database [options]"

    parser = OptionParser(usage)
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=none)")
    parser.add_option("--outpref", dest="outpref", type="string", default="d",
                      help="prefix for output frame or output list with @ mark (default=d)")
    
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()

    distcor(args[0],args[1], inpref=options.inpref, outpref=options.outpref)
 
    
    
