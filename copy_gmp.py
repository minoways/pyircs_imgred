#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
from misc import *
from optparse import OptionParser
import shutil 

def copy_gmp(inlist, inpref='sdfr', outpref='sdf2r'):
    
    # open input list and check if it exists
    ingmp_arr = check_input_file(inlist, inpref, '.gmp')
    if isinstance(ingmp_arr, int):
        return 1

    # check output images 
    outgmp_arr = []
    for i in range(len(ingmp_arr)):
        outgmp_arr.append(ingmp_arr[i].replace(inpref, outpref))

    # check if input images exist and fix bad pixel if requested
    for i in range(len(ingmp_arr)):
        shutil.copy(ingmp_arr[i], outgmp_arr[i])
    
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
        
    copy_gmp(args[0], inpref=options.inpref, outpref=options.outpref)
