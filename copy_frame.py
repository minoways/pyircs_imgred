#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
from misc import *
from optparse import OptionParser

def copy_frame(inlist, inpref='', outpref='r'):
    
    # open input list and check if it exists
    inimg_arr = check_input(inlist, inpref)
    if isinstance(inimg_arr, int):
        return 1

    # check output images 
    outimg_arr = check_outpref(outpref, inimg_arr)
    if isinstance(outimg_arr, int):
        return 1 

    # check if input images exist and fix bad pixel if requested
    iraf.unlearn('imcopy')
    for i in range(len(inimg_arr)):
        iraf.imcopy(inimg_arr[i], outimg_arr[i])
    
    return 0

if __name__=="__main__":

    usage = "usage: %prog input_list [options]"

    parser = OptionParser(usage)
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=null)")
    parser.add_option("--outpref", dest="outpref", type="string", default="r",
                      help="prefix for output frame or output list with @ mark (default=r)")
    
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit()
        
    copy_frame(args[0], inpref=options.inpref, outpref=options.outpref)
