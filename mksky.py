#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
import tempfile 
from zcombine import * 
from misc import *
from optparse import OptionParser

def mksky(inlist, outfits, inpref='', mskpref='none', second=False):
    
    # load IRAF package
    iraf.proto()

    # check input image list 
    inimg_arr = check_input(inlist, inpref)
    if isinstance(inimg_arr, int):
        return 1

    # check output fits
    if os.access(outfits,os.R_OK):
        print >> sys.stderr, 'operation would overwrite exisiting image (%s)' % outfits
        return 1

    # check mask images 
    if mskpref.lower() == 'none' or mskpref == '':
        mask = 0
    else:
        mask = 1
        mskimg_arr = check_inpref(mskpref, inimg_arr)
        if isinstance(mskimg_arr, int):
            return 1

    # check input, output, and mask files and record mask frame into header
    in_subset = []
    if mask == 1:
        mask_subset = []
    iraf.unlearn('hedit')
    for i in range(len(inimg_arr)):
        if mask == 1:
            if mskpref[0] != ':':
                if not second:
                    iraf.hedit(inimg_arr[i], 'objmask', mskimg_arr[i]+'[pl]', add='yes', verify='no')
                else:
                    iraf.hedit(inimg_arr[i], 'objmask', mskimg_arr[i], add='yes', verify='no')

        in_subset.append(inimg_arr[i])
    

    # sky subtraction 
    tmp = tempfile.NamedTemporaryFile(suffix='', prefix='', dir='/tmp')
    tmp_prefix = tmp.name
    tmp.close()

    tmp_subset = '%s_subset_%03d.lst' % (tmp_prefix, int(i))
    # make subset list 
    f = open(tmp_subset, 'w')
    for j in range(len(in_subset)):
        f.write(in_subset[j]+'\n')
    f.close()
    # combine subset 
    if mask == 1:
        ret = zcombine(tmp_subset, outfits, combine='median', reject='sigclip', masktype='!OBJMASK')
    else:
        ret = zcombine(tmp_subset, outfits, combine='median', reject='sigclip', masktype='none')
    if ret != 0:
        remove_temp_all(tmp_prefix)
        return 1

    remove_temp_all(tmp_prefix)
    return 0

if __name__=="__main__":

    usage = "usage: %prog input_list outfits [options]"

    parser = OptionParser(usage)
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=none)")
    parser.add_option("--mskpref", dest="mskpref", type="string", default="none",
                      help="prefix for mask image, image list with \'@\', or header keyword with \':\'")
    parser.add_option("--second", dest="second", action="store_true", default=False,
                      help="Is this second cycle subtraction? (default=False)")
    
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit()
                
    mksky(args[0], args[1], inpref=options.inpref, mskpref=options.mskpref, second=options.second)
  
