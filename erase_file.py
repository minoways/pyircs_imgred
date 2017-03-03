 #!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
from misc import *
from optparse import OptionParser

def erase_file(inlist, inpref='', ext='none'):
    
    # open input list and check if it exists
    try :
        fin = open(inlist)
    except IOError:
        print >> sys.stderr, 'cannot open %s' % inlist
        return 1

    for line in fin:
        if not line.startswith('#'):
            fname = line[:-1]
            if inpref == '':
                inimg = fname
            else:
                inimg = inpref + os.path.basename(fname)
                
            if ext != 'none':
                f,e= os.path.splitext(inimg)
                file_name = f + '.' + ext
                if os.access(file_name, os.R_OK):
                    os.remove(file_name)
            else:
                if os.access(inimg, os.R_OK):
                    os.remove(inimg)
    fin.close()
    
    return 0

if __name__=="__main__":

    usage = "usage: %prog input_list [options]"

    parser = OptionParser(usage)
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix for input frame (default=null)")
    parser.add_option("--ext", dest="ext", type="string", default="none",
                      help="File extension (default=none)")
    
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit()
        
    erase_file(args[0], inpref=options.inpref, ext=options.ext)
