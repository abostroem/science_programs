import pyraf
from pyraf import iraf
from iraf import noao,digiphot,apphot,daofind as daofind
from iraf import noao,digiphot,apphot as apphot
import os
import numpy as np
from astropy.io import fits
import os
from matplotlib import pyplot
import pdb
from matplotlib.backends.backend_pdf import PdfPages
import math
import subprocess
from astropy.io import ascii
from astropy.table import Table
import shutil

#Determine WFC3 PSF

import numpy as np
from scipy.optimize import curve_fit

def fit_func(x, a0, a1, a2, a3, a4, a5):
    '''
    Fits a double gaussian. Used to get the approximate PSF of the WFC3 image
    '''
    z = (x - a1) / a2
    n = (x - a4) / a5
    y = a0 * np.exp(-z**2 / a2) + a3 * np.exp(-n**2 / a5)
    return y

def fit_single_gauss_func(x, a0, a1, a2):
    '''
    Fits a guassian
    '''
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) 
    return y
 
def find_psf(ycenter, xmin):
    '''
    Collapse data along the X axis and fit the profile to get an estimate of the PSF width
    for DAOFIND
    '''
    datax = np.arange(21)+xmin
    datay = img[ycenter, xmin:xmin+21]
    pyplot.plot(datax, datay)
    p0 = [np.max(datay), xmin + 10, 2, np.max(datay), xmin + 10, 2]
    parameters, covariance = curve_fit(fit_func, datax, datay, p0 = p0)
    fitx = np.arange(20000)/1000. + xmin
    fitdata = fit_func(fitx, *parameters)
    pyplot.plot(fitx, fitdata)
    peak = np.max(fitdata)
    halfmax = peak/2.0
    peak_indx = np.where(fitdata == peak)[0][0]
    left_lower_indx = np.where(fitdata[:peak_indx] < halfmax)[0][-1]
    right_lower_indx = np.where(fitdata[peak_indx:] < halfmax)[0][0]
    fwhm = fitx[peak_indx + right_lower_indx] - fitx[left_lower_indx]
    pyplot.plot([fitx[left_lower_indx], fitx[peak_indx + right_lower_indx]], [halfmax, halfmax])
    print 'FWHM (ycenter = %i, xmin = %i) = %f' %(ycenter, xmin, fwhm)

def find_stars_daofind_missed(coord_file, image, ext = 0):
    '''
    Plots the WFC3 image with the points found by the coord_file (from DAOFIND), then the
    user can click on stars that were missed by DAOFIND and the x and y coordinates are appended
    to the DAOFIND coord file
    '''
    with open(coord_file, 'a') as ofile:
        img = fits.getdata(image, ext)
        fig = pyplot.figure()
        ax = fig.add_subplot(1, 1, 1)
        im_plot = ax.imshow(img, interpolation = 'none', cmap = 'bone', vmin = 0, vmax = 50)
        x_coords, y_coords, mag = np.genfromtxt(coord_file, unpack = 'true', usecols = [0, 1, 2])
        ax.plot(x_coords[mag < -4] - 1, y_coords[mag < -4] - 1 , 'o', markerfacecolor = 'none', markeredgecolor = 'lime', markeredgewidth = 2, linestyle = 'none')
        ax.set_ylim(900, 1040)
        add_point = 'y'
        print 'SELECT POINTS FROM THE EDGES'
        raw_input('Zoom in on first missed star')
        while add_point is not 'n':
            coords = pyplot.ginput(n=1)
            x, y = coords[0]
            ax.plot(x, y,'o', markerfacecolor = 'none', markeredgecolor = 'r', markeredgewidth = 2, linestyle = 'none')
            pyplot.draw()
            ofile.write('{}\t{}\n'.format(x+1, y+1))
            add_point = raw_input('Find another point? (y), n')
            ax.plot(x, y, 'o', markerfacecolor = 'none', markeredgecolor = 'lime', markeredgewidth = 2, linestyle = 'none')
            pyplot.draw()
        print 'SELECT POINTS FROM THE CENTER'
        im_plot.set_clim(0, 500)
        pyplot.draw()
        add_point = 'y' 
        while add_point is not 'n':
            coords = pyplot.ginput(n=1)
            x, y = coords[0]
            ax.plot(x, y, 'o',markerfacecolor = 'none', markeredgecolor = 'r', markeredgewidth = 2, linestyle = 'none')
            pyplot.draw()
            ofile.write('{}\t{}\n'.format(x+1, y+1))
            add_point = raw_input('Find another point? (y), n')
            ax.plot(x, y, 'o', markerfacecolor = 'none', markeredgecolor = 'lime', markeredgewidth = 2, linestyle = 'none')
            pyplot.draw()
    
def confirm_daophot(image_file_stis, image_file_wfc3, magh_file_stis, magh_file_wfc3, ext_stis = 0, ext_wfc3 = 0, vmin = 0, vmax = 50):
    '''
    Plots the STIS and WFC3 images side by side and the star locations of stars found by DAOPHOT. 
    '''
    img_stis = fits.getdata(image_file_stis, ext_stis)
    x_coords_stis, y_coords_stis, id_stis, mag = np.genfromtxt(magh_file_stis, usecols = [1, 2, 3, 37], unpack = True, missing = 'INDEF', missing_values = np.nan)
    pp = PdfPages('confirm_daophot_locations_3936.pdf') 
    fig = pyplot.figure(figsize = [25, 15])
    ax_stis = fig.add_subplot(1, 2, 1)
    im_plot_stis = ax_stis.imshow(img_stis, interpolation = 'none', cmap = 'bone', vmin = vmin, vmax = vmax)
    ax_stis.plot(x_coords_stis[(np.isfinite(mag)) & (mag < 19.5)] - 1, y_coords_stis[np.isfinite(mag) & (mag < 19.5)] - 1, 'o', markerfacecolor = 'none', markeredgecolor = 'lime', markeredgewidth = 2, linestyle = 'none')
    for x, y, id in zip(x_coords_stis[(np.isfinite(mag)) & (mag < 19.5)] - 1, y_coords_stis[np.isfinite(mag) & (mag < 19.5)] - 1, id_stis[np.isfinite(mag) & (mag < 19.5)]):
        ax_stis.text(x, y, id, color = 'lime')

    img_wfc3 = fits.getdata(image_file_wfc3, ext_wfc3)
    x_coords_wfc3, y_coords_wfc3, id_wfc3 = np.genfromtxt(magh_file_wfc3, usecols = [1, 2, 3], unpack = True) 
    ax_wfc3 = fig.add_subplot(1, 2, 2)
    im_plot_wfc3 = ax_wfc3.imshow(img_wfc3, interpolation = 'none', cmap = 'bone', vmin = vmin, vmax = vmax)
    ax_wfc3.plot(x_coords_wfc3[np.isfinite(mag) & (mag < 19.5)] - 1, y_coords_wfc3[np.isfinite(mag) & (mag < 19.5)] - 1, 'o', markerfacecolor = 'none', markeredgecolor = 'lime', markeredgewidth = 2, linestyle = 'none')
    for x, y, id in zip(x_coords_wfc3[(np.isfinite(mag)) & (mag < 19.5)] - 1, y_coords_wfc3[np.isfinite(mag) & (mag < 19.5)] - 1, id_wfc3[np.isfinite(mag) & (mag < 19.5)]):
        ax_wfc3.text(x, y, id, color = 'lime')
    pyplot.draw()

    lower_lim = 900
    while lower_lim > 0:
        ax_stis.set_ylim(lower_lim, lower_lim+150)
        ax_wfc3.set_ylim(lower_lim, lower_lim+150)
        pyplot.draw()
        pp.savefig()
        raw_input('Press enter to move down slit')
        lower_lim = lower_lim - 150
    pp.close()
def run_daofind_on_wfc3(image = 'f336w_63.7_cropped.fits', output_find_file = 'daofind_output_wfc3.coo'):
    '''
    Run DAOFIND on the WFC3 image to get coordinates. Produces a file called daofind_output_wfc3.coo. This 
    overwrites a previous file by default
    '''
    
    #Make sure all parameters are default before modifying any of them
    daofind.unlearn()
    apphot.findpars.unlearn()
    apphot.datapars.unlearn()
    #Datapars
    apphot.datapars.scale = 1.0
    apphot.datapars.fwhmpsf = 2.0  #this is a really important parameter
    apphot.datapars.sigma = 0.1
    apphot.datapars.datamin = 0
    apphot.datapars.ccdread = ""
    apphot.datapars.exposure = "EXPTIME"
    #findpars
    apphot.findpars.threshold = 4.0
    apphot.findpars.nsigma = 1.5
    apphot.findpars.theta = 0.0
    apphot.findpars.ratio = 0.5
    apphot.findpars.sharplo = 0.2
    apphot.findpars.sharphi = 10.0
    apphot.findpars.roundlo = -1.0
    apphot.findpars.roundhi = 2.0
    if os.path.exists(os.path.join(os.getcwd(), output_find_file)):
        os.remove(os.path.join(os.getcwd(), output_find_file))
    daofind(image = image, output = output_find_file)

def run_daophot(image, coord_file = 'daofind_output_wfc3.coo', salgorithm = 'median', output_mag_file = None, output_mag_humanreadable = None, apertures = None, maxshift = 1.0):
    '''
    Run DAOPHOT on an image given a coordinate list (coord_file). Apertures, salgorithm, and maxshift are used
    as input parameters to DAOPHOT. This program overwrites an existing file by default
    '''
    if not apertures:
        apertures = '1, 2, 3, 4, 5, 6, 7, 8, 9, 10'
    if not output_mag_file:
        output_mag_file = '{}_{}.mag'.format(os.path.splitext(image)[0], salgorithm)
    if not output_mag_humanreadable:
        output_mag_humanreadable = output_mag_file+'h'
    apphot.centerpars.unlearn()
    apphot.fitskypars.unlearn()
    apphot.fitskypars.salgorithm = salgorithm
    apphot.photpars.apertures = apertures
    apphot.photpars.zmag = 25.0

    if os.path.exists(os.path.join(os.getcwd(), output_mag_file)):
        os.remove(os.path.join(os.getcwd(), output_mag_file))
    apphot.phot(image , coords = coord_file, output = output_mag_file, interactive = 'no', maxshift = maxshift)

    if os.path.exists(os.path.join(os.getcwd(), output_mag_humanreadable)):
        os.remove(os.path.join(os.getcwd(), output_mag_humanreadable))
    apphot.pconvert(textfile = output_mag_file, table = output_mag_humanreadable, fields = 'IMAGE, XCENTER, YCENTER, ID, XINIT, YINIT, RAPERT, SUM, FLUX, MAG, MSKY')


def match_stis_wfc3_stars(image_file_stis, image_file_wfc3, magh_file_stis, magh_file_wfc3, ext_stis = 0, ext_wfc3 = 0):
    '''
    Steps through each object found by DAOPHOT in which the magnitude is not INDEF and is less than 25. The
    Users is asked to verify that the same star was found in both STIS and WFC3 images. If this is the case
    the star is written to an output file. The user is also given the option of changing the contrast
    limits. The output file is called stis_wfc3_coords_mag.dat
    '''
    fig = pyplot.figure(figsize = [15, 10])
    ax_stis = fig.add_subplot(1, 2, 1)
    ax_wfc3 = fig.add_subplot(1, 2, 2)
    img_stis = fits.getdata(image_file_stis, ext_stis)
    x_coords_stis, y_coords_stis, id_stis = np.genfromtxt(magh_file_stis, usecols = [1, 2, 3], unpack = True) 
    img_wfc3 = fits.getdata(image_file_wfc3, ext_wfc3)
    x_coords_wfc3, y_coords_wfc3, id_wfc3, mag = np.genfromtxt(magh_file_wfc3, usecols = [1, 2, 3, 37], unpack = True, missing = 'INDEF', missing_values = np.nan) 
    im_plot_stis = ax_stis.imshow(img_stis, interpolation = 'none', cmap = 'bone', vmin = 0, vmax = 50)
    im_plot_wfc3 = ax_wfc3.imshow(img_wfc3, interpolation = 'none', cmap = 'bone', vmin = 0, vmax = 50)
    indx = mag < 30
    with open('stis_wfc3_coords_mag.dat', 'w') as ofile:
        ofile.write('#ID\tX-STIS\tY-STIS\tX-WFC3\tY-WFC3\tWFC3-Mag-APER= 2\n')
        for ids, idw, xs, xw, ys, yw, m in zip(id_stis[indx], id_wfc3[indx], x_coords_stis[indx], x_coords_wfc3[indx], y_coords_stis[indx], y_coords_wfc3[indx], mag[indx]):
            ax_stis.plot(xs - 1, ys - 1, 'o', markerfacecolor = 'none', markeredgecolor = 'lime', markeredgewidth = 2, linestyle = 'none')
            ax_wfc3.plot(xw - 1, yw - 1, 'o', markerfacecolor = 'none', markeredgecolor = 'lime', markeredgewidth = 2, linestyle = 'none')
            ax_stis.text(xs, ys, ids, color = 'lime', weight = 'bold')
            ax_wfc3.text(xw, yw, idw, color = 'lime', weight = 'bold')
            ax_stis.set_ylim(ys-30, ys+30)
            ax_wfc3.set_ylim(ys-30, ys+30)
            pyplot.draw()
            set_contrast = raw_input('Set new contrast limits for {}? (n), y '.format(ids))
            while set_contrast is 'y':
                vmax = float(raw_input('Enter max limit'))
                im_plot_stis.set_clim(0, vmax)
                im_plot_wfc3.set_clim(0, vmax)
                pyplot.draw()
                set_contrast = raw_input('Set new contrast limits for {}? (y), n '.format(ids))

            write_output = raw_input('Write object {} to output file? (y), n '.format(ids))
            if write_output is not 'n':
                ofile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(ids, xs, ys, xw, yw, m))    
                ax_stis.plot(xs - 1, ys - 1, 'o', markerfacecolor = 'none', markeredgecolor = 'b', markeredgewidth = 2, linestyle = 'none')
                ax_wfc3.plot(xw - 1, yw - 1, 'o', markerfacecolor = 'none', markeredgecolor = 'b', markeredgewidth = 2, linestyle = 'none')
                ax_stis.text(xs, ys, ids, color = 'b', weight = 'bold')
                ax_wfc3.text(xw, yw, idw, color = 'b', weight = 'bold')  
            else:
                ax_stis.plot(xs - 1, ys - 1, 'o', markerfacecolor = 'none', markeredgecolor = 'r', markeredgewidth = 2, linestyle = 'none')
                ax_wfc3.plot(xw - 1, yw - 1, 'o', markerfacecolor = 'none', markeredgecolor = 'r', markeredgewidth = 2, linestyle = 'none')
                ax_stis.text(xs, ys, ids, color = 'r', weight = 'bold')
                ax_wfc3.text(xw, yw, idw, color = 'r', weight = 'bold') 
            ofile.flush()       

def plot_final_stars(image_file_stis, image_file_wfc3, output_file):
    '''
    Plots the final output catalog from make_multispec_output. It steps through the y-axis to give you enough zoom to see the sources.
    The output is saved as final_star_locations_3936.pdf
    '''
    fig = pyplot.figure(figsize = [15, 10])
    ax_stis = fig.add_subplot(1, 2, 1)
    ax_wfc3 = fig.add_subplot(1, 2, 2)
    img_stis = fits.getdata(image_file_stis, 0)
    id, x_coords_stis, y_coords_stis, x_coords_wfc3, y_coords_wfc3, mag = np.genfromtxt(output_file, unpack = True, usecols = [0, 1, 2, 3, 4, 5]) 
    img_wfc3 = fits.getdata(image_file_wfc3, 0)    
    im_plot_stis = ax_stis.imshow(img_stis, interpolation = 'none', cmap = 'bone', vmin = 0, vmax = 200, aspect = 'auto')
    im_plot_wfc3 = ax_wfc3.imshow(img_wfc3, interpolation = 'none', cmap = 'bone', vmin = 0, vmax = 200, aspect = 'auto')
    ax_stis.plot(x_coords_stis -1, y_coords_stis-1, '+', markerfacecolor = 'none', markeredgecolor = 'r', markeredgewidth = 2, linestyle = 'none', markersize = 10) 
    ax_wfc3.plot(x_coords_wfc3-1, y_coords_wfc3-1, '+', markerfacecolor = 'none', markeredgecolor = 'r', markeredgewidth = 2, linestyle = 'none', markersize = 10)
    lower_lim = 900
    pp = PdfPages('final_star_locations_3936.pdf')
    ax_stis.set_title('STIS Psuedo Image')
    ax_wfc3.set_title('WFC3 Cropped Image')
    while lower_lim > 0:
        ax_stis.set_ylim(lower_lim, lower_lim+150)
        ax_wfc3.set_ylim(lower_lim, lower_lim+150)
        pyplot.draw()
        raw_input('Press enter to move down slit')
        pp.savefig()
        lower_lim = lower_lim - 150
    pp.close()
    pyplot.close()

def add_slit(output_file):
    '''
    This function converts the x pixel number to a slit number and overwrites the output file adding the slit number
    '''
    with open(output_file, 'r') as ofile:
        all_lines = ofile.readlines()
    shutil.move(output_file, os.path.splitext(output_file)[0]+'_without_slit'+os.path.splitext(output_file)[1])
    with open(output_file, 'w') as ofile:
        ofile.write('#ID	X-STIS	Y-STIS	X-WFC3	Y-WFC3	WFC3-Mag-APER= 2    SLIT')
        for iline in all_lines[1:]:
            split_line = iline.split()
            slit_num = math.floor(float(split_line[1])/4.) + 1
            ofile.write('{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\n'.format(split_line[0], split_line[1], split_line[2], split_line[3], split_line[4], split_line[5], int(slit_num)))


def test_stis_distortion(output_file):
    '''
    Creates a plot of the distortion of a the STIS CCD by comparing the location of stars in the
    STIS Psuedo image to the location of stars in the WFC3 image. Creates output file
    distortion_3936.pdf
    '''
    pyplot.figure(figsize = [10, 15])
    id, x_coords_stis, y_coords_stis, x_coords_wfc3, y_coords_wfc3, mag = np.genfromtxt(output_file, unpack = True, usecols = [0, 1, 2, 3, 4, 5]) 
    pyplot.plot(x_coords_stis, y_coords_stis, 'bo')
    for xs, ys, xw, yw in zip(x_coords_stis, y_coords_stis, x_coords_wfc3, y_coords_wfc3):
        pyplot.plot([xs, xw]*50, [ys, yw]*50, 'r')
    pyplot.xlabel('X pixel')
    pyplot.ylabel('Y pixel')
    pyplot.title('CCD distortion of G430M/3936 median y offset = {} pixels'.format(np.median(abs(y_coords_stis - y_coords_wfc3))))
    pyplot.legend(['STIS Psuedo Image', 'WFC3 Image to STIS Psuedo Image * 50'], loc = 'best')
    raw_input('Enter to continue')
    pyplot.savefig('distortion_3936.pdf')
    pyplot.close()

def check_for_drift_during_visit(output_file):
    '''
    Creates a plot to look for drift between exposures (slit positions) by looking at the
    change in y location with slit number. A linear trend is fit. Figure saved as drift_3936.pdf
    '''
    id, x_coords_stis, y_coords_stis, x_coords_wfc3, y_coords_wfc3, mag, slit = np.genfromtxt(output_file, unpack = True)
    pyplot.figure()
    pyplot.plot(slit, y_coords_wfc3 - y_coords_stis, 'o')
    fit_coeff_all = np.polyfit(slit, y_coords_wfc3 - y_coords_stis, 1)
    pyplot.plot(np.arange(18), np.polyval(fit_coeff_all, np.arange(18)))
    fit_coeff_se = np.polyfit(slit[slit < 10], y_coords_wfc3[slit< 10] - y_coords_stis[slit < 10], 1)
    pyplot.plot(np.arange(10), np.polyval(fit_coeff_se, np.arange(10)))
    fit_coeff_nw = np.polyfit(slit[slit >= 10], y_coords_wfc3[slit >= 10] - y_coords_stis[slit >= 10], 1)
    pyplot.plot(np.arange(10, 18, 1), np.polyval(fit_coeff_nw, np.arange(10, 18, 1)))
    pyplot.axvline(10.5)
    pyplot.legend(['Data', 'y = {}x + {}'.format(fit_coeff_all[0],  fit_coeff_all[1]), 'y = {}x + {}'.format(fit_coeff_se[0],  fit_coeff_se[1]), 'y = {}x + {}'.format(fit_coeff_nw[0],  fit_coeff_nw[1]), '6 month split'])
    pyplot.xlabel('X Pixel (slit)')
    pyplot.ylabel('Difference in Y location between STIS and WFC3')
    pyplot.draw()
    raw_input('Press Enter to finish')
    pyplot.savefig('drift_3936.pdf')
    apply_drift = raw_input('Apply drift correction? ')
    if apply_drift == 'y':
        shutil.move(output_file, os.path.splitext(output_file)[0]+'_uncorrected'+os.path.splitext(output_file)[1])
        y_coords_stis[slit<10] = y_coords_stis[slit<10] + fit_coeff_se[0]*slit[slit<10] + fit_coeff_se[1]
        y_coords_stis[slit>=10] = y_coords_stis[slit>=10] + fit_coeff_nw[0]*slit[slit>=10] + fit_coeff_nw[1]
        table_data = Table([id, x_coords_stis, y_coords_stis, x_coords_wfc3, y_coords_wfc3, mag, slit], names = ['ID', 'X-STIS', 'Y-STIS', 'X-WFC3', 'Y-WFC3', 'WFC3-Mag-APER= 2', 'SLIT'])
        ascii.write(table_data, output_file, format = 'commented_header')

def make_offset_histograms(output_file):
    '''
    Creates a histogram of the difference in star location in terms of x and in terms of y.
    A guassian is fit to the histogram to estimate the mean and standard deviation
    '''
    id, x_coords_stis, y_coords_stis, x_coords_wfc3, y_coords_wfc3, mag = np.genfromtxt(output_file, unpack = True, usecols = [0, 1, 2, 3, 4, 5]) 
    pyplot.figure()
    n, bins, patches = pyplot.hist(x_coords_wfc3 - x_coords_stis, bins = np.arange(-10, 10, 0.5))
    parameters, covariance = curve_fit(fit_single_gauss_func, bins[:-1], n)  
    fitdata = fit_single_gauss_func(np.arange(-10, 10, 0.1), *parameters)
    pyplot.plot(np.arange(-10, 10, 0.1), fitdata)
    pyplot.title('Delta X (WFC3 - STIS), mean = {:2.4}, stdev = {:2.4}'.format(parameters[1], parameters[2]))
    pyplot.xlabel('Delta X')
    pyplot.savefig('delta_x_3936.pdf')
    pyplot.close()

    pyplot.figure()
    n, bins, patches = pyplot.hist(y_coords_wfc3 - y_coords_stis, bins = np.arange(-5, 5, 0.05))
    parameters, covariance = curve_fit(fit_single_gauss_func, bins[:-1], n)  
    fitdata = fit_single_gauss_func(np.arange(-5, 5, 0.05), *parameters)
    pyplot.plot(np.arange(-5, 5, 0.05), fitdata)
    pyplot.title('Delta Y (WFC3 - STIS) mean = {:2.4}, stdev = {:2.4}'.format(parameters[1], parameters[2]))
    pyplot.xlabel('Delta Y')
    pyplot.savefig('delta_y_3936.pdf')
    pyplot.close()


def make_input_match_catalog(output_file, wfc3_img, wfpc2_img):
    id_wfc3, x_coords_wfc3, y_coords_wfc3, mag_wfc3, slit = np.genfromtxt(output_file, unpack = True, usecols = [0, 3, 4, 5, 6])
    id_wfpc2, x_coords_wfpc2, y_coords_wfpc2, mag_wfpc2 = np.genfromtxt('/Users/bostroem/science/multispec/multi_spec_files/hunter_ApJ1995_488_179.dat', unpack = True, usecols = [0, 1, 2, 3])
    indx1 = np.argsort(mag_wfpc2)
    indx2 = [(x_coords_wfpc2[indx1] > 300) & (x_coords_wfpc2[indx1] < 430) & (mag_wfpc2[indx1] < 90.)]
    indx = indx1[indx2]
    fig = pyplot.figure(figsize = [15, 15])
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    im1 = ax1.imshow(fits.getdata(wfc3_img, 0), cmap = 'bone', vmin = 0, vmax = 50)
    im2 = ax2.imshow(fits.getdata(wfpc2_img, 1), cmap = 'bone', vmin = 0, vmax = 50)
    ax1.plot(x_coords_wfc3-1, y_coords_wfc3-1, 'r.')
    ax2.plot(x_coords_wfpc2[indx]-1, y_coords_wfpc2[indx]-1, 'r.')
    ax2.set_xlim(300, 430)
    for x, y, id in zip(x_coords_wfc3, y_coords_wfc3, id_wfc3):
        ax1.text(x, y, id, color = 'c')
    for x, y, id in zip(x_coords_wfpc2[indx], y_coords_wfpc2[indx], id_wfpc2[indx]):
        ax2.text(x, y, id, color = 'c')
    ax1.set_ylim(700, 900)
    ax2.set_ylim(600, 800)
    pyplot.draw()
    with open('geomap_input.dat', 'w') as ofile:
        print 'Enter n to exit'
        wfc3_id = 'y'
        while wfc3_id is not 'n':
            try:
                wfc3_id = float(raw_input('Enter the WFC3 star ID '))
                wfpc2_id = float(raw_input('Enter the WFPC2 star ID '))
                wfc3_indx = [id_wfc3 == wfc3_id]
                wfpc2_indx = [id_wfpc2 == wfpc2_id]
                ofile.write('{}\t{}\t{}\t{}\n'.format(x_coords_wfc3[wfc3_indx][0], y_coords_wfc3[wfc3_indx][0], x_coords_wfpc2[wfpc2_indx][0], y_coords_wfpc2[wfpc2_indx][0]))
            except (IndexError, ValueError) as e:
                print e.message
                wfc3_id = 'n'
        pyplot.close()


def get_nearest_obj(x_coords_wfpc2, y_coords_wfpc2, x_wfpc2, y_wfpc2, mag_wfpc2):
    dist = np.sqrt((x_coords_wfpc2 - x_wfpc2)**2 + (y_coords_wfpc2 - y_wfpc2)**2)
    max_dist = 2
    nearby_indx = np.where(dist < 2)[0]
    if np.min(dist) > 2:
        return None
    while len(nearby_indx) == 0:
        max_dist +=1
        nearby_indx = np.where(dist < max_dist)[0]
        pdb.set_trace()
    brightest_indx = np.argmin(mag_wfpc2[nearby_indx])
    return nearby_indx[brightest_indx]


def get_hunter_id_numbers(output_file, wfc3_img, wfpc2_img):
    #make_input_match_catalog(output_file, wfc3_img, wfpc2_img)
    import iraf
    from pyraf import iraf
    from iraf import images,immatch,geomap as geomap
    from iraf import images,immatch,geoxytran as geoxytran
    id_wfc3, x_coords_wfc3, y_coords_wfc3, mag_wfc3, slit = np.genfromtxt(output_file, unpack = True, usecols = [0, 3, 4, 5, 6])
    make_input_match_catalog(output_file, wfc3_img, wfpc2_img)
    geomap('geomap_input.dat', 'geomap_output.txt', 0, 70, 0, 1024, interactive = False)
    data = Table([x_coords_wfc3, y_coords_wfc3], names = ['x', 'y'])
    ascii.write(data, 'wfc3_coord_list.tab', format = 'no_header')
    if os.path.exists(os.path.join(os.getcwd(), 'wfpc2_coord_from_geoxytran.dat')):
        os.remove(os.path.join(os.getcwd(), 'wfpc2_coord_from_geoxytran.dat'))
    geoxytran(input = 'wfc3_coord_list.tab', output = 'wfpc2_coord_from_geoxytran.dat', database = 'geomap_output.txt', transforms = 'geomap_input.dat')
    x_coords_wfpc2_guess, y_coords_wfpc2_guess = np.genfromtxt('wfpc2_coord_from_geoxytran.dat', unpack = True)
    fig = pyplot.figure()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    im1 = ax1.imshow(fits.getdata(wfc3_img, 0), cmap = 'bone', vmin = 0, vmax = 50)
    im2 = ax2.imshow(fits.getdata(wfpc2_img, 1), cmap = 'bone', vmin = 0, vmax = 50)
    ax1.plot(x_coords_wfc3-1, y_coords_wfc3-1, 'r.')
    ax2.plot(x_coords_wfpc2_guess-1, y_coords_wfpc2_guess-1, 'r.')
    ax2.set_xlim(300, 430)
    for x, y, id in zip(x_coords_wfc3, y_coords_wfc3, id_wfc3):
        ax1.text(x, y, id, color = 'c')
    for x, y, id in zip(x_coords_wfpc2_guess, y_coords_wfpc2_guess, id_wfc3):
        ax2.text(x, y, id, color = 'c')
    hunter_id = []
    hunter_mag = []
    hunter_x_coords = []
    hunter_y_coords = []
    id_wfpc2, x_coords_wfpc2, y_coords_wfpc2, mag_wfpc2 = np.genfromtxt('/Users/bostroem/science/multispec/multi_spec_files/hunter_ApJ1995_488_179.dat', unpack = True, usecols = [0, 1, 2, 3])
    for i, x_wfpc2, y_wfpc2 in zip(range(len(x_coords_wfc3)), x_coords_wfpc2_guess, y_coords_wfpc2_guess):
        if (y_wfpc2 < 40) | (y_wfpc2 > 800): #Outside the range of the WFPC2 image
            hunter_id.append('J{}'.format(id_wfc3[i]))
            hunter_mag.append(mag_wfc3[i])
            hunter_x_coords.append(0)
            hunter_y_coords.append(0)
        else:
            indx = get_nearest_obj(x_coords_wfpc2, y_coords_wfpc2, x_wfpc2, y_wfpc2, mag_wfpc2)
            if indx:
                hunter_id.append('H{}'.format(id_wfpc2[indx]))
                hunter_mag.append(mag_wfpc2[indx])
                hunter_x_coords.append(x_coords_wfpc2[indx])
                hunter_y_coords.append(y_coords_wfpc2[indx])
            else:
                hunter_id.append('J{}'.format(id_wfc3[i]))
                hunter_mag.append(mag_wfc3[i])
                hunter_x_coords.append(0)
                hunter_y_coords.append(0)
    hunter_id = np.array(hunter_id)
    hunter_mag = np.array(hunter_mag)
    hunter_x_coords = np.array(hunter_x_coords)
    hunter_y_coords = np.array(hunter_y_coords)
    output_table = Table([hunter_id, x_coords_wfc3, y_coords_wfc3, hunter_x_coords, hunter_y_coords, hunter_mag, slit], names = ['id', 'wfc3_x', 'wfc3_y', 'wfpc2_x', 'wfpc2_y', 'mag', 'slit'])
    ascii.write(output_table, 'final_list_w_hunter_id.dat', format = 'fixed_width')


def plot_final_star_match(image_file_wfc3, image_file_wfpc2, coord_file):
    '''
    Plots the final output catalog from make_multispec_output. It steps through the y-axis to give you enough zoom to see the sources.
    The output is saved as final_star_locations_3936.pdf
    '''
    fig = pyplot.figure(figsize = [15, 10])
    ax_wfc3 = fig.add_subplot(1, 2, 1)
    ax_wfpc2 = fig.add_subplot(1, 2, 2)
    img_wfc3 = fits.getdata(image_file_wfc3, 0)
    table_data = ascii.read(coord_file) 
    x_coords_wfc3 = table_data['wfc3_x'].data
    y_coords_wfc3 = table_data['wfc3_y'].data
    x_coords_wfpc2 = table_data['wfpc2_x'].data
    y_coords_wfpc2 = table_data['wfpc2_y'].data
    id = table_data['id'].data
    img_wfpc2 = fits.getdata(image_file_wfpc2, 1)    
    im_plot_wfc3 = ax_wfc3.imshow(img_wfc3, interpolation = 'none', cmap = 'bone', vmin = 0, vmax = 200, aspect = 'auto')
    im_plot_wfpc2 = ax_wfpc2.imshow(img_wfpc2, interpolation = 'none', cmap = 'bone', vmin = 0, vmax = 40, aspect = 'auto')
    ax_wfc3.plot(x_coords_wfc3 -1, y_coords_wfc3-1, '+', markerfacecolor = 'none', markeredgecolor = 'r', markeredgewidth = 2, linestyle = 'none', markersize = 10) 
    ax_wfpc2.plot(x_coords_wfpc2-1, y_coords_wfpc2-1, '+', markerfacecolor = 'none', markeredgecolor = 'r', markeredgewidth = 2, linestyle = 'none', markersize = 10)
    for i, x, y in zip(id, x_coords_wfc3, y_coords_wfc3):
        ax_wfc3.text(x, y, i, color = 'r')
    pyplot.draw()
    for i, x, y in zip(id, x_coords_wfpc2, y_coords_wfpc2):
        ax_wfpc2.text(x, y, i, color = 'r')
    pyplot.draw()
    lower_lim = 900
    ax_wfpc2.set_xlim(335, 425)
    pp = PdfPages('star_match_wfc3_wfpc2_locations_3936.pdf')
    ax_wfc3.set_title('WFC3 Cropped Image')
    ax_wfpc2.set_title('WFPC2 Image')
    while lower_lim > 0:
        print lower_lim
        ax_wfc3.set_ylim(lower_lim, lower_lim+150)
        ax_wfpc2.set_ylim((min(705,lower_lim-160))*1.1, (lower_lim-160+150)*1.1)
        pyplot.draw()
        raw_input('Press enter to move down slit')
        pp.savefig()
        lower_lim = lower_lim - 150
    pp.close()
    pyplot.close()

def generate_multispec_input(coord_file, slit_num):
    table_data = ascii.read(coord_file)
    indx = np.where(table_data['slit'].data == slit_num)[0]
    e = np.ones((len(indx),))*0.400
    x = np.ones((len(indx),))*511.5
    ra = np.zeros((len(indx),))
    dec = np.zeros((len(indx),))
    sed = ['kurucz_30000_ms_z0.fits']*len(indx)
    slit_id = table_data['id'].data[indx]
    slit_y = table_data['wfc3_y'].data[indx]-1.  #Change from DAOFIND indexing (center of pixel at 1) to IDL indexing (center of pixel at 0)
    slit_mag = table_data['mag'].data[indx]
    data = Table([slit_id, x, slit_y, ra, dec, sed, e, slit_mag], names = ['ID', 'x', 'y', 'RA', 'DEC', 'sed', 'e(4405-5495)', 'wfpc2_f336w'])
    if slit_num < 10:
        ascii.write(data, 'slit0{}_phot.dat'.format(slit_num), format = 'tab')
    else:
        ascii.write(data, 'slit{}_phot.dat'.format(slit_num), format = 'tab')
    
def print_to_screen(what_to_print):
    print '#-----------------------------'
    print what_to_print
    print '#-----------------------------'

if __name__ == "__main__":
    os.chdir('/Users/bostroem/science/multispec/ccd_multispec')
    salgorithm = 'centroid'
    #Determine WFC3 PSF
    img = fits.getdata('f336w_63.7_cropped.fits', 0)
    #-----------------------------
    #Find PSF
    #-----------------------------
    #print '#-----------------------------'
    #print 'Analyzing PSF'
    #sprint '#-----------------------------'
    #find_psf(654, 0)
    #find_psf(924, 16)
    #find_psf(315, 32)
    #find_psf(6, 24)
    #find_psf(172, 32)

    #-----------------------------
    #Create Coordinate List
    #-----------------------------
    print_to_screen('run DAOFIND on WFC3 Image')
    #run_daofind_on_wfc3()
    print_to_screen('Find stars DAOFIND missed')
    #find_stars_daofind_missed('daofind_output_wfc3.coo', 'f336w_63.7_cropped.fits' , ext = 0)
    print_to_screen('run DAOPHOT on STIS and WFC3 images')
    #Get photometry and STIS image locations
    #run_daophot('f336w_63.7_cropped.fits', salgorithm = salgorithm)
    #run_daophot('long_slit_img_3936.fits', maxshift = 4.0, salgorithm = salgorithm)

    #-----------------------------
    #Check photometry results
    #-----------------------------
    #print_to_screen("Inspect where DAOPHOT found stars, duplicates will be removed later")
    #confirm_daophot('long_slit_img_3936.fits', 'f336w_63.7_cropped.fits', 'long_slit_img_3936_{}.magh'.format(salgorithm), 'f336w_63.7_cropped_{}.magh'.format(salgorithm), vmin = 0, vmax = 800)
    #print 'with a different contrast'
    #confirm_daophot('long_slit_img_3936.fits', 'f336w_63.7_cropped.fits', 'long_slit_img_3936_{}.magh'.format(salgorithm), 'f336w_63.7_cropped_{}.magh'.format(salgorithm))

    #-----------------------------
    #Create Multispec Output file
    #-----------------------------
    #print_to_screen('ID good matches between WFC3 and STIS stars')
    #match_stis_wfc3_stars('long_slit_img_3936.fits', 'f336w_63.7_cropped.fits', 'long_slit_img_3936_{}.magh'.format(salgorithm), 'f336w_63.7_cropped_{}.magh'.format(salgorithm))

    #-----------------------------
    #Check final locations
    #-----------------------------
    #print_to_screen('Plot final selection of stars')
    #plot_final_stars('long_slit_img_3936.fits', 'f336w_63.7_cropped.fits', 'stis_wfc3_coords_mag.dat')

    #-----------------------------
    #Add in the slit number to the Multispec file
    #-----------------------------
    #print_to_screen('Add slit number to output file')
    #add_slit('stis_wfc3_coords_mag.dat')

    #-----------------------------
    #Check for drifts from exposure to exposure
    #-----------------------------
    #print_to_screen('Correct for drift within a visit')
    #check_for_drift_during_visit('stis_wfc3_coords_mag.dat')

    #-----------------------------
    #Create a histogram of the offset in X and Y
    #-----------------------------
    #print_to_screen('Look for offsets between the images in x and y')
    #make_offset_histograms('stis_wfc3_coords_mag.dat')

    #-----------------------------
    #Test distortion
    #-----------------------------
    #print_to_screen('Check for geometric distortion')
    #test_stis_distortion('stis_wfc3_coords_mag.dat')
    #-----------------------------
    #Make multispec input files
    #-----------------------------

    print_to_screen('Match WFC3 stars to WFPC2 Hunter stars')
    get_hunter_id_numbers('stis_wfc3_coords_mag.dat', 'f336w_63.7_cropped.fits', '/Users/bostroem/science/multispec/multi_spec_files/r136_f555w_wfpc2_images/u25y0105t_c0m.fits')
    plot_final_star_match('f336w_63.7_cropped.fits', '/Users/bostroem/science/multispec/multi_spec_files/r136_f555w_wfpc2_images/u25y0105t_c0m.fits', 'final_list_w_hunter_id.dat')
    print_to_screen('Make an _phot.dat file for each slit')
    for i in range(1, 17, 1):
        generate_multispec_input('final_list_w_hunter_id.dat', i)


