#!/usr/bin/env python

import os,sys
import pyfits
#from pyds9 import *
from ds9 import ds9 
import stsci.numdisplay as nd
from optparse import OptionParser

def frcheck(infile):

    #<--- count number of lines in the list and backup the orignal file --->

    # open list
    try :
        f=open(infile)
    except IOError:
        print >> sys.stderr, 'cannot open %s' % infile
        return 1

    nline = 0
    orgfile = infile+".org"
    f_org = open(orgfile,'w')
    for line in f:
        f_org.write(line)
        if not line.startswith('#'):
            nline += 1
    # close list 
    f.close()
    f_org.close()

    #<--- input the file names into array --->#
    # open list 
    f = open(infile)
    fname=[[] for i in range(nline)]
    frchk=[[] for i in range(nline)]
    i = 0
    for line in f:
        if not line.startswith('#'):
            fname[i] = line[:-1]
            if not os.access(fname[i], os.R_OK):
                print >> sys.stderr, "cannot read %s" % (fname[i])
                f.close()
                return 1
            i += 1
    # close list
    f.close()

    #<--- check frame --->
    for i in range(nline):
        im1 = pyfits.open(fname[i])
        coadds = float(im1[0].header['COADDS'])
        ndr = float(im1[0].header['NDR'])
        base1,ext = os.path.splitext(os.path.basename(fname[i]))
        if i != nline - 1:
            im2 = pyfits.open(fname[i+1])
            base2,ext = os.path.splitext(os.path.basename(fname[i+1]))
        else :
            im2 = pyfits.open(fname[i-1])
            base2,ext = os.path.splitext(os.path.basename(fname[i-1]))
        
        name = '%s-%s' % (base1,base2)
        sub = im1[0].data - im2[0].data
        
        z1_tmp = -1000.0 * float(coadds) * float(ndr)
        z2_tmp = 1000.0 * float(coadds) * float(ndr)

        zs = nd.zscale.zscale(sub,nsamples=1000, contrast=0.25)
        #print 'test: %s %s' % (str(zs[0]),str(zs[1])) 
        if str(zs[0]) != 'nan':
            z1_tmp = zs[0]
        if str(zs[1]) != 'nan':
            z2_tmp = zs[1]
        
        nd.display(sub, bufname='imt1024', name=name, z1=z1_tmp, z2=z2_tmp)

        d.set('zoom to fit')
        check = ''
        while check != 'OK' and check != 'NG':
            check = raw_input('%s - %s (ok/ng): ' %(base1, base2))
            if check.lower() == '' or check.lower() == 'ok':
                check = 'OK'
            elif check.lower() == 'ng':
                check = 'NG'
            else :
                print 'Error: unknown answer (%s)' % (check)
        print check
        if check == 'OK':
            frchk[i] = 0
        else: 
            frchk[i] = 1

        # close fits files
        im1.close()
        im2.close()
        
        
    #<--- output bad removed list --->#
    # open list
    try:
        f=open(orgfile)
    except IOError:
        print >> sys.stderr, 'cannot read original file %s' % (orgfile)
    f_out = open(infile, 'w')
    i = 0
    for line in f:
        if not line.startswith('#'):
            if i == 0:
                if frchk[i] == 0:
                    f_out.write('%s\n' % fname[i])
                else :
                    f_out.write('#%s\n'% fname[i])
                i += 1
            else:
                if frchk[i] == 0 or frchk[i-1] == 0:
                    f_out.write('%s\n' % fname[i])
                else :
                    f_out.write('#%s\n'% fname[i])
                i += 1
        else:
            f_out.write(line)

    # close list 
    f.close()
    f_out.close()

    return 0


if __name__=="__main__":
    argvs=sys.argv
    argc=len(argvs)

    if argc < 2:
        print 'Usage: %s file_lst' %(os.path.basename(sys.argv[0]))
        sys.exit()
        
    d = ds9()
    #d = DS9()
    frcheck(sys.argv[1])
