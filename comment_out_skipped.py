#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
from misc import *
from optparse import OptionParser
import shutil
import tempfile 

def comment_out_skipped(inlist, inpref='', sspref='sdfr'):
    
    # open list
    try :
        f=open(inlist)
    except IOError:
        print >> sys.stderr, 'cannot open %s' % inlist
        return 1

    # prepare for temporary file 
    tmp = tempfile.NamedTemporaryFile(suffix='', prefix='', dir='/tmp')
    tmp_prefix = tmp.name
    tmp.close()

    rawfile = tmp_prefix + '_rawfile.dat'
    orgfile = inlist+".skiporg"
    f_org = open(orgfile,'w')
    f_out = open(rawfile,'w')
    for line in f:
        if not line.startswith('#'):
            fname, ext = os.path.splitext(os.path.basename(line[:-1]))
            gmp = sspref + fname + '.gmp'

            try:
                fgmp = open(gmp)
            except:
                print >> sys.stderr, 'cannot read geomap file: %s' % gmp
                remove_temp_all(tmp_prefix)
                return 1

            nobj = 0
            for gline in fgmp:
                if not gline.startswith('#'):
                    nobj += 1
            fgmp.close()

            if nobj > 0:
                f_out.write(line)
            else:
                f_out.write('#'+line)
        else:
            f_out.write(line)

        f_org.write(line)
    # close list 
    f.close()
    f_org.close()
    f_out.close()
    
    # copy 
    shutil.copy(rawfile, inlist)
    # remove temporary files
    remove_temp_all(tmp_prefix)

    return 0

if __name__=="__main__":

    usage = "usage: %prog input_list [options]"

    parser = OptionParser(usage)
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=none)")
    parser.add_option("--sspref", dest="sspref", type="string", default="sdfr",
                      help="prefix for ss frame or output list with @ mark (default=sdfr)")
    
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit()
        
    comment_out_skipped(args[0], inpref=options.inpref, sspref=options.sspref)
