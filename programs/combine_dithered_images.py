import pyfits
import numpy as np
import os
import glob
import c_correlate
import scipy
from scipy import signal
from scipy.signal import medfilt 
from scipy import stats
from scipy.stats import tstd

import pdb
from matplotlib import pyplot
from optparse import OptionParser

from pyraf import iraf
from iraf import stsdas,hst_calib,stis,ocrreject as ocrreject

#!!!!!!!!!!!!!!!!!!!!!!
#CORRECTNESS OF DQ AND ERR ARRAYS HAVE NOT BEEN CONFIRMED ALTHOUGH THEY SHOULD FOLLOW THE SAME
#RULES AS THE IMAGE EXTENSION
#!!!!!!!!!!!!!!!!!!!!!!


def make_declination_dict(flist):
    '''
    Create a dictionary with declination as keys and filenames as lists to each key
    This is to match 2 dithered images together

    Create a dictionary with declination as keys and filenames as lists to each key
    This is to match 2 dithered images together

    Called from: Main
    Calls to: Nothing

    Input:
        flist: list of files to make dictionary out of
    Outputs:
        Nothing
    Returns:
        dec_dict: a dictionary with declination as keys and filenames as lists to each key
    '''
    dec_dict = {} 
    for ifile in flist:
        targ_dec = pyfits.getval(ifile, 'dec_targ', 0)
        if dec_dict.has_key(targ_dec):
            dec_dict[targ_dec].append(ifile)
        else:
            dec_dict[targ_dec] = [ifile]
    return dec_dict

def combine_dithered_images(dec_dict, targ_dec, use_hdr_offset):
    '''
    The 12465 data is dithered along the slit direction. This program finds the shifts
    of the dithered images by cross-correlating a middle column and combines them.

    Input:
        dec_dict: dictionary with declination as keys and a list of filenames of data taken at 
            the that declination as values
        targ_dec: decclination of target
        offset: dither offset in pixels
    Outputs:
        None
    Returns:
        new_img: combined image array
    '''

    #Read in 2 images
    ofile1 = pyfits.open(dec_dict[targ_dec][0])
    ofile2 = pyfits.open(dec_dict[targ_dec][1])
    hdr0 = ofile1[0].header 
    hdr1 = ofile1[1].header
    hdr2 = ofile1[2].header
    hdr3 = ofile1[3].header
    img1 = ofile1[1].data
    err1 = ofile1[2].data
    dq1 = ofile1[3].data
    img2 = ofile2[1].data
    err2 = ofile2[2].data
    dq2 = ofile2[3].data
    hdr02 = ofile2[0].header
    postarg1 = hdr0['postarg2']
    postarg2 = hdr02['postarg2']
    detector = hdr0['detector']
    if detector == 'FUV-MAMA':
        plate_scale = 0.0246
    if detector == 'CCD':
        plate_scale = 0.05078
    max_lag = postarg2*plate_scale - postarg1*plate_scale

    cenwave = hdr0['cenwave']
    exptime2 = pyfits.getval(dec_dict[targ_dec][1], 'texptime', 0)
    hdr0['texptime'] = hdr0['texptime'] + exptime2  #exptime is only in the Primary HDU
    if use_hdr_offset is False:
        #Get middle columns of each image
        col1 = np.sum(img1[:, 350:650], axis = 1)
        col2 = np.sum(img2[:, 350:650], axis = 1)
        #Calculate offset
        lag, corr_mat = c_correlate.c_corr(col1, col2)
        max_indx = np.argmax(corr_mat)
        print 'Found offset of %i' %(lag[max_indx])
        max_lag = lag[max_indx]
        if max_lag == 0:
            lag = np.append(lag[:max_indx], lag[max_indx+1:])
            corr_mat = np.append(corr_mat[:max_indx], corr_mat[max_indx+1:])
            max_indx = np.argmax(corr_mat)
            print 'Found new offset of %i' %(lag[max_indx])
    else:
        postarg1 = hdr0['postarg2']
        postarg2 = hdr02['postarg2']
        detector = hdr0['detector']
        if detector == 'FUV-MAMA':
            plate_scale = 0.0246
        if detector == 'CCD':
            plate_scale = 0.05078
        max_lag = int(postarg2/plate_scale - postarg1/plate_scale)
        print 'Using offset of %i found from postarg2 header keyword' %(max_lag)

    if max_lag > 0:
        start1 = 0
        end1 = np.shape(img1)[0]
        start2 = max_lag
        end2 = np.shape(img1)[0] + max_lag
        if not use_hdr_offset:
            pdb.set_trace()
        create_new_combined_arrays(img1, img2, err1, err2, dq1, dq2, start1, end1, start2, end2, ofile1[0].header, ofile1[1].header, ofile1[2].header, ofile1[3].header, ofile2[0].header, ofile2[1].header, ofile2[2].header, ofile2[3].header)

    else:
        start1 = 0
        end1 = np.shape(img1)[0]
        start2 = np.abs(max_lag)
        end2 = np.shape(img1)[0]+np.abs(max_lag)
        if not use_hdr_offset:
            pdb.set_trace()
        create_new_combined_arrays(img2, img1, err2, err1, dq2, dq1, start1, end1, start2, end2, ofile2[0].header, ofile2[1].header, ofile2[2].header, ofile2[3].header, ofile1[0].header, ofile1[1].header, ofile1[2].header, ofile1[3].header)

    #Write fits file
    #make_combined_fits(hdr0, hdr1, hdr2, hdr3, new_img, new_err, new_dq, cenwave)
 
   

def id_cr(orig_img, kernel_size = 9, thresh = 300):
    img = orig_img.copy()
    filt_img = img.copy()
    l1 = pyplot.imshow(img, interpolation = 'nearest', cmap = 'bone', vmin = 0, vmax = 1000)
    for irow in range(np.shape(img)[0]):
        filt_img[irow, :] = medfilt(filt_img[irow, :], kernel_size = kernel_size)
        stdev = tstd(img[irow, :])
        cr_pix = np.where(abs((img[irow, :] - filt_img[irow, :])) > stdev)
        img[irow, cr_pix[0]] = -999
    #What should I do with the error array
    x = np.where(img == -999)
    #Uncomment this block of code to see 
    #print np.shape(x)
    #pyplot.plot(x[1], x[0], 'r.')
    #pdb.set_trace()
    #pyplot.close()

    return img

def cr_reject(new_img1, new_img2):
    #Identify cosmic rays
    cr_id_img1 = id_cr(new_img1)
    cr_id_img2 = id_cr(new_img2)
    #If a good pixel exists in one image, replace the bad pixel with the good pixel
    img1_replace_indx = np.where((cr_id_img1 == -999) & (cr_id_img2 != -999))
    cr_id_img1[img1_replace_indx] = cr_id_img2[img1_replace_indx]
    img2_replace_indx = np.where((cr_id_img2 == -999) & (cr_id_img1 != -999))
    #Set pixels with are bad in both images to 0
    cr_id_img1[cr_id_img1 == -999] = 0.0
    cr_id_img2[cr_id_img2 == -999] = 0.0
    new_img = cr_id_img1 + cr_id_img2
    return new_img

def create_new_combined_arrays(top_img, bottom_img, top_err, bottom_err, top_dq, bottom_dq, start1, end1, start2, end2, hdr10, hdr11, hdr12, hdr13, hdr20, hdr21, hdr22, hdr23):
        new_img1 = np.zeros((end2, end1))
        new_img2 = np.zeros((end2, end1))
        new_img1[start1:end1, :] = new_img1[start1:end1, :] + bottom_img
        new_img2[start2:end2, :] = new_img2[start2:end2, :]+  top_img        

        new_err1 = np.zeros((end2, end1))
        new_err2 = np.zeros((end2, end1))
        new_err1[start1:end1, :] = new_err1[start1:end1, :] + bottom_err
        new_err2[start2:end2, :] = new_err2[start2:end2, :] + top_err
        new_err = np.sqrt(0.25*new_err1**2 + 0.25*new_err2**2)

        new_dq1 = np.zeros((end2, end1))
        new_dq2 = np.zeros((end2, end1))
        new_dq1[start1:end1, :] = new_dq1[start1:end1, :] + bottom_dq
        new_dq2[start2:end2, :] = new_dq2[start2:end2, :] + top_dq
        new_dq1 = np.int_(new_dq1)
        new_dq2 = np.int_(new_dq2)
        new_dq = np.bitwise_or(new_dq1, new_dq2)
        
        if hdr10['detector'] == 'CCD':
            make_combined_fits(hdr10, hdr11, hdr12, hdr13, new_img1, new_err1, new_dq1, suffix = 'ocr_tmp1.fits')
            make_combined_fits(hdr20, hdr21, hdr22, hdr23, new_img2, new_err2, new_dq2, suffix = 'ocr_tmp2.fits')
            if os.path.exists('%s_%i_combined_img.fits' %(hdr10['targname'][5:], hdr10['cenwave'])):
                os.remove('%s_%i_combined_img.fits' %(hdr10['targname'][5:], hdr10['cenwave']))
            ocrreject(input = ','.join(glob.glob('%s_ocr_tmp*.fits' %(hdr10['targname'][5:]))), output = '%s_%i_combined_img.fits' %(hdr10['targname'][5:], hdr10['cenwave']))
        else:
            new_img = (new_img1+new_img2)/2.0
            make_combined_fits(hdr10, hdr11, hdr12, hdr13, new_img, new_err, new_dq, cenwave = hdr10['cenwave'])
        


def make_combined_fits(hdr0, hdr1, hdr2, hdr3, new_img, new_err, new_dq, cenwave = None, suffix = 'combined_img.fits'):
    '''
    write combined image to file
    Called from:
        combine_dithered_images
    Calls to:
        Nothing
    Inputs:
        hdr: fits header object
        img: image array
        cenwave: central wavelength for filename
    Outputs:
        TARGNAME_(CENWAVE)_combined_img.fits files of the combined dithered images
    Returns:
        Nothing

    '''
    hdu0 = pyfits.PrimaryHDU()
    hdu0.header = hdr0
    print '#############################%s######################' %(hdr0['targname'])
    hdu1 = pyfits.ImageHDU(data=new_img, name='SCI', header = hdr1)
    hdu2 = pyfits.ImageHDU(data=new_err, name='ERR', header = hdr2)
    hdu3 = pyfits.ImageHDU(data = np.int16(new_dq), name = 'DQ', header = hdr3)
    hdulist = pyfits.HDUList([hdu0, hdu1, hdu2, hdu3])
    if cenwave:
        hdulist.writeto('%s_%i_%s' %(hdr0['targname'][5:], hdr0['cenwave'], suffix), clobber = True)
    else:
        hdulist.writeto('%s_%s' %(hdr0['targname'][5:], suffix), clobber = True)



if __name__ == "__main__":

    os.environ['oref'] = '/Users/bostroem/science/oref/'
    parser = OptionParser()
    parser.add_option('--UseHeader', dest = 'use_hdr_offset', action = 'store_true', help = 'Set the dithered offset in pixels', default = False)
    (options, args) = parser.parse_args()
    '''
    #for the FUV data
    idir = '/user/bostroem/science/12465_otfr20120425/mama/'

    ##idir = '/Users/bostroem/science/12465_otfr20120425/mama/'
    os.chdir(idir)
    flist = glob.glob('obrc04???_flt.fits')+glob.glob('obrc05???_flt.fits')
    dec_dict = make_declination_dict(flist)
    for targ_dec in dec_dict.keys():
        combine_dithered_images(dec_dict, targ_dec, options.use_hdr_offset)
    pdb.set_trace()
    '''
    #idir = '/Users/bostroem/science/12465_otfr20121109/ccd/'
    #idir = '/user/bostroem/science/12465_otfr20130503/ccd/'
    idir = '/Users/bostroem/science/2014_dc_aas/'
    os.chdir(idir)
    iraf.chdir(idir)
    

    '''    
    flist = glob.glob('obrc06*_flt.fits')#+glob.glob('ob???????_flt.fits')
    dec_dict = make_declination_dict(flist)
    for targ_dec in dec_dict.keys():
        combine_dithered_images(dec_dict, targ_dec, options.use_hdr_offset)
    '''    
    flist = glob.glob('obrc01*_flc_flt.fits')#+glob.glob('ob???????_flc_flt.fits')
    dec_dict = make_declination_dict(flist)
    for targ_dec in dec_dict.keys():
        combine_dithered_images(dec_dict, targ_dec, options.use_hdr_offset)
    '''
    flist = glob.glob('obrc07*_flt.fits')#+glob.glob('ob???????_flt.fits')
    dec_dict = make_declination_dict(flist)
    for targ_dec in dec_dict.keys():
        combine_dithered_images(dec_dict, targ_dec, options.use_hdr_offset)
    '''
    flist = glob.glob('obrc02*_flc_flt.fits')#+glob.glob('ob???????_flc_flt.fits')
    dec_dict = make_declination_dict(flist)
    for targ_dec in dec_dict.keys():
        combine_dithered_images(dec_dict, targ_dec, options.use_hdr_offset)
    '''
    flist = glob.glob('obrc08*_flt.fits')#+glob.glob('ob???????_flt.fits')
    dec_dict = make_declination_dict(flist)
    for targ_dec in dec_dict.keys():
        combine_dithered_images(dec_dict, targ_dec, options.use_hdr_offset)
    '''
    flist = glob.glob('obrc03*_flc_flt.fits')#+glob.glob('ob???????_flc_flt.fits')
    dec_dict = make_declination_dict(flist)
    for targ_dec in dec_dict.keys():
        combine_dithered_images(dec_dict, targ_dec, options.use_hdr_offset)
    '''
    flist = glob.glob('obrc09*_flt.fits')#+glob.glob('ob???????_flt.fits')
    dec_dict = make_declination_dict(flist)
    for targ_dec in dec_dict.keys():
        combine_dithered_images(dec_dict, targ_dec, options.use_hdr_offset)

    flist = glob.glob('obzk01*_flt.fits')#+glob.glob('ob???????_flt.fits')
    dec_dict = make_declination_dict(flist)
    for targ_dec in dec_dict.keys():
        combine_dithered_images(dec_dict, targ_dec, options.use_hdr_offset)
    '''


    #obrc01, obrc07: 3936
    #obrc02, obrc08: 4451
    #obrc03, obrc09: 4706
    #obzk01: 4194
    #obzk02: 4194

