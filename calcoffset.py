#!/usr/bin/env python
import os,sys
import pyfits
from pyraf import iraf
import tempfile 
from zcombine import * 
from misc import *
from ds9 import ds9 
from optparse import OptionParser

def calcoffset(inlist, inpref='', trace=True, review=False):

    # open input list and check if it exists
    inimg_arr = check_input(inlist, inpref)
    if isinstance(inimg_arr,int):
        return 1 

    # check output geomap file
    gmp_arr = []
    for i in range(len(inimg_arr)):
        # geomap input file name 
        fname,ext = os.path.splitext(inimg_arr[i])
        gmp_arr.append(fname + '.gmp')

        if os.access(gmp_arr[i], os.R_OK):
            print >> sys.stderr, 'operation would overwrite existing file (%s)' % gmp_arr[i]
            return 1


    # open ds9 
    d = ds9()
    d.set('regions delete all')

    # unlearn IRAF commands
    iraf.reset(stdimage='imt1024')
    iraf.unlearn('display')
    iraf.unlearn('rimexam')
    iraf.unlearn('imexam')
    iraf.unlearn('hedit')
    
    # review rimexam parameters
    if review:
        iraf.epar('rimexam')

    # prefix for temporary file 
    tmp = tempfile.NamedTemporaryFile(suffix='', prefix='', dir='/tmp')
    tmp_prefix = tmp.name
    tmp.close()

    # save imexam parameters into temporary files
    tmp_par = tmp_prefix + '.par'
    iraf.rimexam.saveParList(filename=tmp_par)
    
    # region file for the first image
    tmp0_reg = tmp_prefix+'_0.reg'

    # measure position 
    i = 0
    k = 0
    xarr = []
    yarr = []
    peak_arr = []
    fwhm_arr = []
    gsx0 = 0.0
    gsy0 = 0.0
    nobj = 1
    trace0 = trace
    for i in range(len(inimg_arr)):
        # initialize trace parameter 
        trace = trace0

        # get dithering position from header 
        im = pyfits.open(inimg_arr[i], mode='update')
        try:
            ao_mode = im[0].header['D_MODE']
            ao_loop = im[0].header['D_LOOP']
        except KeyError:
            trace = False

        # read guide star coordinates
        gsx = 0.0
        gsy = 0.0
        if ao_loop.lower() == 'on':
            if ao_mode.lower().find('lgs') == -1:
                try:
                    gsx = float(im[0].header['D_AU1GSX'])
                    gsy = float(im[0].header['D_AU1GSY'])
                except KeyError:
                    trace = False
            else:
                try:
                    gsx = float(im[0].header['D_AU2GSX'])
                    gsy = float(im[0].header['D_AU2GSY'])
                except KeyError:
                    trace = False
        else:
            trace = False

        if k == 0:
            gsx0 = gsx
            gsy0 = gsy

        # load rimexam parameters
        iraf.rimexam.setParList(ParList=tmp_par)

        # display image
        d.set('regions delete all')
        iraf.display(inimg_arr[i], 1)

        # automatic object pickup
        calc_ng = 1
        if k >0 and trace == True:
            xobj = []
            yobj = []
            peak = []
            fwhm = []

            # guess object position 
            tmp_reg = tmp_prefix+'.reg'
            if os.access(tmp_reg, os.R_OK):
                os.remove(tmp_reg)
            freg = open(tmp_reg, 'w')
            for ii in range(nobj):
                xg = xref[ii] + (gsx - gsx0)
                yg = yref[ii] + (gsy - gsy0)            
                freg.write('image; circle %.3f %.3f 10 # color=blue text={%d}\n' % (xg, yg, (ii+1)))
            freg.close()
            d.set('regions load %s' % tmp_reg)

            # check object position guess
            OutOfRange = 0
            for ii in range(nobj):
                xg = xref[ii] + (gsx - gsx0)
                yg = yref[ii] + (gsy - gsy0)
                if xg <= 0 or xg >= 1024 or yg <= 0 or yg >= 1024:
                    OutOfRange = 1
            
            # pickup objects
            if OutOfRange == 0:
                if os.access(tmp_reg, os.R_OK):
                    os.remove(tmp_reg)
                freg = open(tmp_reg, 'w')
                for ii in range(nobj):
                    xg = xref[ii] + (gsx - gsx0)
                    yg = yref[ii] + (gsy - gsy0)            

                    tmp_coo = tmp_prefix+'_coo.dat'
                    if os.access(tmp_coo, os.R_OK):
                        os.remove(tmp_coo)
                    fcoo = open(tmp_coo, 'w')
                    fcoo.write('%.3f %.3f a\n' % (xg, yg))
                    fcoo.close()
                    imexam_ng = 1
                    auto_skip = 0
                    while imexam_ng == 1:
                        try:
                            ret = iraf.imexam(inimg_arr[i],1,imagecu=tmp_coo, use_dis='no', Stdout=1)
                            imexam_ng = 0
                            auto_skip = 0
                        except:
                            print '\nIRAF imexam failed for object No. %d' % (ii+1) 
                            print 'Try again by changing rimexam parameters'
                            print 'Hit return to enter rimexam parameter setting window'
                            print 'Type \'q\' to skip this object.'
                    
                            check = ''
                            while check.lower() != 'q' and check.lower() != 'rimexam':
                                check = raw_input('Hit return or type \'q\':')
                                if check.lower() == '' or check.lower() == 'rimexam':
                                    check = 'rimexam'
                                    print 'Push Save&Quit button to quit from the parameter setting window'
                                    iraf.epar('rimexam')
                                    imexam_ng = 1
                                elif check.lower() == 'q':
                                    check = 'q'
                                    auto_skip = 1
                                    imexam_ng = 0
                                else :
                                    print 'Error: unknown answer (%s)' % (check)
                
                    os.remove(tmp_coo)

                    if auto_skip == 0:
                        # display result
                        for j in range(len(ret)):
                            print ret[j]
                        print '\n'

                        # parse results
                        param1 = ret[len(ret)-2].split()
                        param2 = ret[len(ret)-1].split()
                        if len(param1) == 4 and len(param2) == 11:
                            if isfloat(param1[0]):
                                xobj.append(float(param1[0]))
                            if isfloat(param1[0]):
                                yobj.append(float(param1[1]))
                            if isfloat(param2[4]):
                                peak.append(float(param2[4]))
                            else:
                                peak.append(-9999.0)
                            if isfloat(param2[9]):
                                fwhm.append(float(param2[9]))
                            else:
                                fwhm.append(-9999.0)
                            freg.write('image; circle %.3f %.3f 5 # color=red\n' % (xobj[ii], yobj[ii]))
                
                        else :
                            xobj.append(-9999.0)
                            yobj.append(-9999.0)
                            peak.append(-9999.0)
                            fwhm.append(-9999.0)
                    else:
                        xobj.append(-9999.0)
                        yobj.append(-9999.0)
                        peak.append(-9999.0)
                        fwhm.append(-9999.0)

                freg.close()
                d.set('regions load %s' % tmp_reg)
                os.remove(tmp_reg)

                check = raw_input('Is this position okay? (yes/no/skip)')
                if check.lower() == '' or check.lower() == 'yes':
                    calc_ng = 0

                    # update reference points
                    nskip = 0
                    for ii in range(len(xobj)):
                        if xobj[ii] < 0:
                            nskip += 1
                    if nobj == len(xobj) and nskip == 0:
                        gsx0 = gsx
                        gsy0 = gsy 
                        xref = xobj
                        yref = yobj 
                elif check.lower() == 'skip':
                    calc_ng = 0
                    skip_image = 1
                    xobj = []
                    yobj = []
                    peak = []
                    fwhm = []
                    for ii in range(nobj):
                        xobj.append(-9999.0)
                        yobj.append(-9999.0)
                        peak.append(-9999.0)
                        fwhm.append(-9999.0)
                else :
                    calc_ng = 1
            else:
                print ""
                print "Warning:"
                print "Estimated position is out of range"
                print "Pick up object manually\n"
                calc_ng = 1
                
        # calculate position and fwhm
        if calc_ng == 1:

            skip_image = 0
            imexam_ng = 1
            calc_ok = 0
            while calc_ok != 1:
                # display image
                d.set('regions delete all')
                iraf.display(inimg_arr[i], 1)
                # show position of the objects in the first image
                if k != 0:
                    if not trace:
                        d.set('regions load %s' % tmp0_reg)

                    xobj = []
                    yobj = []
                    peak = []
                    fwhm = []
                    skip_obj = 0
                    n_skip_obj = 0
                    for nn in range(nobj):
                        imexam_ng  = 1
                        skip_obj = 0

                        if trace:
                            xg = xref[nn] + (gsx - gsx0)
                            yg = yref[nn] + (gsy - gsy0)            

                            tmp_reg = tmp_prefix+'.reg'
                            if os.access(tmp_reg, os.R_OK):
                                os.remove(tmp_reg)
                            freg = open(tmp_reg, 'w')
                            freg.write('image; circle %.3f %.3f 20 # color=blue text={%d}\n' % (xg, yg, (nn+1)))
                            freg.close()
                            d.set('regions delete all')
                            d.set('regions load %s' % tmp_reg)

                        while imexam_ng == 1:
                            print '\n'
                            print '##### Pickup object No. %d / %d #####' % (nn+1, nobj)
                            print '\n'
                            print '\nSelect Object --> type a , Quit --> type q on the ds9 image'
                            print 'To skip this image type q on the ds9 image and then this image will not be used'
                            print 'Caution: do not type any key except a or q'
                            try :
                                ret = iraf.imexam(inimg_arr[i],1,Stdout=1)
                                if len(ret) == 0:
                                    skip_obj = 1
                                while len(ret) < 4 and skip_obj == 0:
                                    print '\nSelect Object --> type a , Quit --> type q on the ds9 image'
                                    print 'To skip this image type q on the ds9 image and then this image will not be used'
                                    print 'Caution: do not type any key except a or q'
                                    ret = iraf.imexam(inimg_arr[i],1,Stdout=1)
                                    if len(ret) == 0:
                                        skip_obj = 1
                                imexam_ng = 0
                            except:
                                print '\nIRAF imexam failed' 
                                print 'Try again by changing rimexam parameters'
                                print 'Hit return to enter rimexam parameter setting window'
                                print 'To skip this image type \'q\' and then this image will not be used.'
                            
                                check = ''
                                while check.lower() != 'q' and check.lower() != 'rimexam':
                                    check = raw_input('Hit return or type \'q\':')
                                    if check.lower() == '' or check.lower() == 'rimexam':
                                        check = 'rimexam'
                                        print 'Push Save\&Quit to quit from the parameter setting window'
                                        iraf.epar('rimexam')
                                        imexam_ng = 1
                                    elif check.lower() == 'q':
                                        check = 'q'
                                        skip_obj = 1
                                        imexam_ng = 0
                                    else :
                                        print 'Error: unknown answer (%s)' % (check)
                        
                        if skip_obj == 0:
                            # display result
                            for j in range(len(ret)):
                                print ret[j]
                            print '\n'
            
                            # parse results
                            for ii in range(len(ret)):
                                if not ret[ii].startswith('#') and len(ret[ii].split()) > 2:
                                    param = ret[ii].split()
                                    if len(param) == 4:
                                        if isfloat(param[0]):
                                            xobj_tmp = float(param[0])
                                        if isfloat(param[1]):
                                            yobj_tmp = float(param[1])
                                    elif len(param) == 11:
                                        if isfloat(param[4]):
                                            peak_tmp = float(param[4])
                                        else:
                                            peak_tmp = -9999.0
                                        if isfloat(param[9]):
                                            fwhm_tmp = float(param[9])
                                        else:
                                            fwhm_tmp = -9999.0
                                    else:
                                        print >> sys.stderr, 'failed to pick up object in %s' % inimg
                                        remove_temp_all(tmp_prefix)
                                        return 1
                            xobj.append(xobj_tmp)
                            yobj.append(yobj_tmp)
                            peak.append(peak_tmp)
                            fwhm.append(fwhm_tmp)
                        else:
                            xobj.append(-9999.0)
                            yobj.append(-9999.0)
                            peak.append(-9999.0)
                            fwhm.append(-9999.0)
                            n_skip_obj += 1
                        
                        # check number of skipped objects 
                        if n_skip_obj == nobj:
                            skip_image = 1
                            calc_ok = 1
                        else:
                            skip_image = 0
                else:
                    xobj = []
                    yobj = []
                    peak = []
                    fwhm = []
                    while imexam_ng == 1:
                        print '\nSelect Object --> type a , Quit --> type q on the ds9 image'
                        print 'To skip this image type q on the ds9 image and then this image will not be used'
                        print 'Caution: do not type any key except a or q'
                        try :
                            ret = iraf.imexam(inimg_arr[i],1,Stdout=1)
                            if len(ret) == 0:
                                skip_image = 1
                            while len(ret) < 4 and skip_image == 0:
                                print '\nSelect Object --> type a , Quit --> type q on the ds9 image'
                                print 'To skip this image type q on the ds9 image and then this image will not be used'
                                print 'Caution: do not type any key except a or q'
                                ret = iraf.imexam(inimg_arr[i],1,Stdout=1)
                                if len(ret) == 0:
                                    skip_image = 1
                            imexam_ng = 0
                        except:
                            print '\nIRAF imexam failed' 
                            print 'Try again by changing rimexam parameters'
                            print 'Hit return to enter rimexam parameter setting window'
                            print 'To skip this image type \'q\' and then this image will not be used.'
                            
                            check = ''
                            while check.lower() != 'q' and check.lower() != 'rimexam':
                                check = raw_input('Hit return or type \'q\':')
                                if check.lower() == '' or check.lower() == 'rimexam':
                                    check = 'rimexam'
                                    print 'Push Save\&Quit to quit from the parameter setting window'
                                    iraf.epar('rimexam')
                                    imexam_ng = 1
                                elif check.lower() == 'q':
                                    check = 'q'
                                    skip_image = 1
                                    imexam_ng = 0
                                else :
                                    print 'Error: unknown answer (%s)' % (check)

                    if skip_image == 0:
                        # display result
                        for j in range(len(ret)):
                            print ret[j]
                        print '\n'
            
                        # parse results
                        for ii in range(len(ret)):
                            if not ret[ii].startswith('#') and len(ret[ii].split()) > 2:
                                param = ret[ii].split()
                                if len(param) == 4:
                                    if isfloat(param[0]):
                                        xobj.append(float(param[0]))
                                    if isfloat(param[1]):
                                        yobj.append(float(param[1]))
                                elif len(param) == 11:
                                    if isfloat(param[4]):
                                        peak.append(float(param[4]))
                                    else:
                                        peak.append(-9999.0)
                                    if isfloat(param[9]):
                                        fwhm.append(float(param[9]))
                                    else:
                                        fwhm.append(-9999.0)
                                else:
                                    print >> sys.stderr, 'failed to pick up object in %s' % inimg
                                    remove_temp_all(tmp_prefix)
                                    return 1

                        # check consistency 
                        nobj = len(xobj)
                        if nobj != len(yobj) or nobj != len(peak) or nobj != len(fwhm):
                            print >> sys.stderr, 'Number of the recorded objects is inconsistent'
                            remove_temp_all(tmp_prefix)
                            return 1

                        xref = xobj
                        yref = yobj
                        
                        # save dummy values for the skipped frames at the beginning 
                        if i != 0:
                            for ii in range(i):
                                xarr[ii] = []
                                yarr[ii] = []
                                for jj in range(nobj):
                                    xarr[ii].append(-9999.0)
                                    yarr[ii].append(-9999.0)
                                    peak_arr[ii].append(-9999.0)
                                    fwhm_arr[ii].append(-9999.0)

                        # save position of the objects in the first image
                        if os.access(tmp0_reg, os.R_OK):
                            os.remove(tmp0_reg)
                        freg0 = open(tmp0_reg, 'w')
                        for ii in range(nobj):
                            freg0.write('image; point(%.3f,%.3f) # point=x text={%d}\n' % (xobj[ii], yobj[ii], (ii+1)))
                        freg0.close()
                    else:
                        xobj.append(-9999.0)
                        yobj.append(-9999.0)
                        peak.append(-9999.0)
                        fwhm.append(-9999.0)
                        calc_ok = 1

                if skip_image == 0:
                    # show results on ds9
                    tmp_reg = tmp_prefix+'.reg'
                    if os.access(tmp_reg, os.R_OK):
                        os.remove(tmp_reg)
                    freg = open(tmp_reg, 'w')
                    if len(xobj) > 0:
                        for ii in range(nobj):
                            if xobj[ii] >= 0:
                                freg.write('image; circle %.3f %.3f 5 # color=red text={%d}\n' % (xobj[ii], yobj[ii], (ii+1)))
                    freg.close()
                    d.set('regions load %s' % tmp_reg)
                    os.remove(tmp_reg)

                    print 'Is this position okay? '
                    print '<options>'
                    print ' yes or return : accept this postion'
                    print ' no or n       : measure position again'
                    print ' r or rimexam  : change rimexam parameters'
                    print ' q or quit     : skip this image'
                    check = raw_input('')
                    if check.lower() == '' or check.lower() == 'yes':
                        calc_ok = 1
                        xref = xobj
                        yref = yobj
                        gsx0 = gsx
                        gsy0 = gsy
                        
                    elif check.lower() == 'r' or check.lower() == 'rimexam':
                        print 'Push Save\&Quit to quit from the parameter setting window\n'
                        iraf.epar('rimexam')
                        calc_ok = 0
                        imexam_ng = 1
                    elif check.lower() == 'q' or check.lower() == 'quit':
                        print 'Skip image (%s)' % inimg_arr[i]
                        xobj = []
                        yobj = []
                        peak = []
                        fwhm = []
                        if k == 0:
                            xobj.append(-9999.0)
                            yobj.append(-9999.0)
                            peak.append(-9999.0)
                            fwhm.append(-9999.0)
                        else:
                            for ii in range(nobj):
                                xobj.append(-9999.0)
                                yobj.append(-9999.0)
                                peak.append(-9999.0)
                                fwhm.append(-9999.0)
                        calc_ok = 1
                    else : 
                        imexam_ng = 1
                        calc_ok = 0

        # save position into array 
        xarr.append(xobj)
        yarr.append(yobj)
        peak_arr.append(peak)
        fwhm_arr.append(fwhm)

        # close image handler 
        im.close()

        # increment counter
        if skip_image == 0:
            k += 1
        
    # remove all temporary files
    remove_temp_all(tmp_prefix)
    
    # save measured coordinates into geomap input file and fits header 
    for i in range(len(inimg_arr)):

        # open image handler 
        im = pyfits.open(inimg_arr[i], mode='update')

        # open file handler 
        fgmp = open(gmp_arr[i], 'w')

        for j in range(nobj):
            # record imexam results into header
            key = 'xc%d' % (j+1)
            im[0].header.update(key,xarr[i][j])
            key = 'yc%d' % (j+1)
            im[0].header.update(key,yarr[i][j])
            key = 'peak%d' % (j+1)
            im[0].header.update(key,peak_arr[i][j])
            key = 'fwhm%d' % (j+1)
            im[0].header.update(key,fwhm_arr[i][j])

            if xarr[i][j] >= 0:
                fgmp.write('%.3f %.3f %.3f %.3f\n' % (xref[j],yref[j],xarr[i][j],yarr[i][j]))
            else:
                fgmp.write('#%.3f %.3f %.3f %.3f\n' % (xref[j],yref[j],xarr[i][j],yarr[i][j]))

        # close image handler
        im.close()

        # close file handler 
        fgmp.close()

    return 0


if __name__=="__main__":

    usage = "usage: %prog input_list [options]"

    parser = OptionParser(usage)
    parser.add_option("--inpref", dest="inpref", type="string", default="",
                      help="prefix of input frame (default=none)")
    parser.add_option("--trace", dest="trace", action="store_true", default=False,
                      help="automatic trace of object position (default=True)")
    parser.add_option("--review", dest="review", action="store_true", default=False,
                      help="review imexam parameters (default=False)")
    
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit()
    
    calcoffset(args[0], inpref=options.inpref, trace=options.trace, review=options.review)
  
