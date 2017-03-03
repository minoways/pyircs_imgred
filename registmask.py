#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
from misc import *
from optparse import OptionParser

def registmask(inlist, inpref='', mskpref='mc', key='OBJMASK'):
    
    # open input list and check if it exists
    inimg_arr = check_input(inlist, inpref)
    if isinstance(inimg_arr, int):
        return 1

    # check output images 
    mskimg_arr = check_inpref(mskpref, inimg_arr)
    if isinstance(mskimg_arr, int):
        return 1 

    # check if input images exist and fix bad pixel if requested
    iraf.unlearn('hedit')
    for i in range(len(inimg_arr)):
        iraf.hedit(inimg_arr[i], key, mskimg_arr[i], add='yes', verify='no')
    
    return 0

if __name__=="__main__":

    usage = "usage: %prog input_list [options]"

    parser = OptionParser(usage)
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=none)")
    parser.add_option("--mskpref", dest="mskpref", type="string", default="mc",
                      help="prefix for mask frame or output list with @ mark (default=mc)")
    parser.add_option("--key", dest="key", type="string", default="OBJMASK",
                      help="header keyword for mask frame (default=OBJMASK)")
    
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit()
        
    registmask(args[0], inpref=options.inpref, mskpref=options.mskpref, key=options.key)
