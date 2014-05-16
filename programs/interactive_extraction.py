import pyfits
import numpy as np
import sys
import math
from matplotlib import pyplot
import pdb
import os
import shutil
from optparse import OptionParser
from scipy.interpolate import LSQUnivariateSpline


#import pyraf
#from pyraf import iraf
#from iraf import stsdas,hst_calib,stis,x1d

import stistools
from stistools.x1d import x1d

def collapse_spectrum(img, num_cols):
    '''Collapse 2D spectrum in the dispersion direction to create a cross dispersion profile'''
    collapsed_img = np.sum(img[:, int(512 - num_cols/2.0):int(512 + num_cols/2.0)], axis = 1)/num_cols
    return collapsed_img

def get_background_locations(input, fig1, ax1, pix_num, collapsed_img, cenwave, first_time = True):
    '''Interactively select the regions to be fit as background '''
    finished_flag = 'n'
    background_collapsed = np.empty((0,))
    background_pix = np.empty((0,))
    #no background regions defined
    if first_time is True:
        if os.path.exists(input.filename.replace('.fits', '_c%i_background_regions.txt' %(cenwave))):
            keep_file_flag = raw_input('Background file already exists: Overwrite file (w) or Append to file (a)') #you should never enter this loop if mode = passive
        else:
            keep_file_flag = 'w'
        ofile = open(input.filename.replace('.fits','_c%i_background_regions.txt' %(cenwave)), keep_file_flag)
    #Background regions have been defined - adding or removing some of them
    if first_time is False:
        remove_flag = 'y'
        #Get initial indx values from background file
        start_indx, end_indx  = np.genfromtxt(input.filename.replace('.fits','_c%i_background_regions.txt' %(cenwave)), unpack = True)
        sort_indx = np.argsort(start_indx)
        start_indx = start_indx[sort_indx]
        end_indx = end_indx[sort_indx]
        #REMOVE BACKGROUND REGIONS
        while remove_flag != 'n':
            l2 = [] #list for line objects (this is so we can make them invisible later)
            l3 = [] #list of text objects (this is so we can make them invisible later)
            #Draw background regions as currently defined
            for i, sindx, eindx in zip(range(len(start_indx)), start_indx, end_indx):
                if sindx < eindx:
                    temp_indx = np.where((pix_num <= eindx) & (pix_num >= sindx))[0]
                else:  #in case you defined a background region from right to left
                    temp_indx = np.where((pix_num <= sindx) & (pix_num >= eindx))[0]
                temp_pix = pix_num[temp_indx]
                temp_bkg = collapsed_img[temp_indx]
                l2.append(pyplot.plot([pix_num[temp_indx[0]], pix_num[temp_indx[-1] + 1]], [collapsed_img[temp_indx[0]], collapsed_img[temp_indx[-1]+1]], 'mo--', markersize = 5))
                l3.append(pyplot.text(np.mean(temp_pix), np.max(temp_bkg), str(i)))
            remove_flag = interact_w_user(input.mode, message = 'Would you like to remove any background regions? (y), n ', default = 'n')
            if remove_flag == 'n':
                break
            try:
                remove_indx = int(interact_w_user(input.mode, message = 'Enter the number of the region you wish to remove: ', default = 1)) #you should never enter this loop if mode = passive
            except ValueError:
                print 'You must enter something that can be converted into an integer'
                remove_indx = int(interact_w_user(input.mode, message = 'Enter the number of the region you wish to remove: ', default = 1)) #you should never enter this loop if mode = passive
            #remove indx from start_indx and end_indx
            start_indx = np.append(start_indx[0:remove_indx], start_indx[remove_indx+1:])
            end_indx = np.append(end_indx[0:remove_indx], end_indx[remove_indx+1:])
            #delete lines from plot
            for iline, itxt in zip(l2, l3):
                iline[0].set_visible(False)
                itxt.set_visible(False)
        ofile = open(input.filename.replace('.fits','_c%i_background_regions.txt' %(cenwave)), 'w')
        for sindx, eindx in zip(start_indx, end_indx):
            temp_indx = np.where((pix_num <= eindx) & (pix_num >= sindx))[0]
            background_collapsed = np.append(background_collapsed, collapsed_img[temp_indx[0]:temp_indx[-1]+1])
            background_pix = np.append(background_pix, pix_num[temp_indx[0]:temp_indx[-1]+1])
            ofile.write('%i\t %i \n' %(int(sindx),int(eindx)))
        finished_flag = interact_w_user(input.mode, message = 'Finished entering points? (n), y ', default = 'y')

    #ADD BACKGROUND REGIONS
    while finished_flag != 'y':
        print 'Select background points for cross-dispersion background subtraction'
        temp_region = fig1.ginput(n = 2, timeout = -1)
        back_reg = [temp_region[0][0], temp_region[1][0]]
        if first_time is True:
            l1 = pyplot.plot(back_reg, np.interp(np.array(back_reg), pix_num, collapsed_img), '|-', color = '#990000',  markersize = 5)
        else:
            l1 = pyplot.plot(back_reg, np.interp(np.array(back_reg), pix_num, collapsed_img), '|-', markersize = 5)
        pyplot.draw()
        keep_flag = interact_w_user(input.mode, message = 'Keep points? (y), n ', default = 'y')
        if keep_flag != 'n':
            ofile.write('%i\t %i \n' %(int(math.floor(temp_region[0][0])),int(math.ceil(temp_region[1][0]))))
            background_collapsed = np.append(background_collapsed, collapsed_img[int(math.floor(temp_region[0][0])):int(math.ceil(temp_region[1][0]))+1])
            background_pix = np.append(background_pix, pix_num[int(math.floor(temp_region[0][0])):int(math.ceil(temp_region[1][0]))+1]) #This only works if you select L to R, I should fix that
            #pdb.set_trace()
        else:
            l1[0].set_visible(False)
            pyplot.draw()
        finished_flag = interact_w_user(input.mode, message = 'Finished entering points? (n), y', default = 'y')
    ofile.close()
    sort_indx = np.argsort(background_pix)
    background_pix = background_pix[sort_indx]
    background_collapsed = background_collapsed[sort_indx]
    #pdb.set_trace()
    try:  #remove any residual background line labeling
        for iline, itxt in zip(l2, l3):
            iline[0].set_visible(False)
            itxt.set_visible(False)
    except:
        pass
    try:
        l1[0].set_visible(False)
    except:
        pass
    return fig1, ax1, background_pix, background_collapsed

def fit_background(input, fig1, ax1, background_pix, background_collapsed, pix_num):
    '''Interactively fit background regions with a spline '''
    change_deg_flag = 'y'
    leg_lines = []
    leg_text = []
    for i in np.arange(5)+1:
        fit = spline_fit(background_pix, background_collapsed, pix_num, ax1, order = i)
        l = ax1.plot(pix_num, fit)
        leg_lines.append(l[0])
        leg_text.append('spline, order %i' %(i))
    leg = ax1.legend(leg_lines, leg_text)
    fig1.canvas.draw()
    change_deg_flag = interact_w_user(input.mode, message = 'Select which order spline you would like to use for your fit (1, 2, 3, 4, 5), -or- redefine background, r ', default = 3)
    if change_deg_flag != 'r':
        while int(change_deg_flag) not in [1, 2, 3, 4, 5]:
            change_deg_flag = raw_input('Enter a valid spline degree: ')  #you should never enter this loop if mode = passive
        fit = spline_fit(background_pix, background_collapsed, pix_num, ax1, order = int(change_deg_flag))
    else:
        for l in leg_lines:
            l.set_visible(False)
        return None, None
    return ax1, fit

def spline_fit(background_pix, background_collapsed, pix_num, ax1, order = 3):
    '''Fit background w/ Spline'''
    sub_values = background_pix[1:] - background_pix[:-1]
    stop_indx = np.append(np.where(sub_values > 1.0)[0], -1)
    start_indx = np.append(0, np.where(sub_values > 1.0)[0] + 1)

    #Place nodes at the begining and ending of each background region
    nodes = np.append(background_pix[start_indx[1:]], background_pix[stop_indx[:-1]])
    nodes = list(set(nodes))  #make sure each node is unique
    nodes.sort()
    y = LSQUnivariateSpline(background_pix, background_collapsed, nodes, k = order)
    ###For debugging, plots the node locations and the points IDed as background regions###
    #ax1.plot(background_pix, background_collapsed, '.')
    #for n in nodes: ax1.axvline(n, color = 'r', linestyle = ':')
    return y(pix_num)



def subtract_2D_cross_disp_background(img, fit):
    '''Create a 2D image from the 1D fit to the background in the cross dispersion direction and subtract it from the image '''
    two_D_background = np.rot90(np.rot90(np.rot90(np.tile(fit, (1024, 1)))))
    return img - two_D_background

def write_temp_subtracted_file(input, subtracted_img):
    '''Write the subtracted image to a temporary file for calstis to extract from. This file is erased at the end of execution '''
    shutil.copyfile(input.filename, input.filename.replace('.fits', 'sub.fits'))
    ofile = pyfits.open(input.filename.replace('.fits', 'sub.fits'), mode = 'update')
    ofile[1].data = subtracted_img
    ofile.flush()
    ofile.close()

def select_extraction_location(fig1, c, region_to_mark = 'spectrum'):  #you should never enter this function if mode = passive
    '''Interactvely select the center of an extraction or background region '''
    interact_w_user(input.mode, message = 'Zoom in on spectrum, hit Enter when finished')
    print 'Select %s location' %(region_to_mark)
    keep_flag = 'n'
    while keep_flag == 'n':
        coords = fig1.ginput(n = 1, timeout = -1)
        ln = pyplot.plot([coords[0][0]], [coords[0][1]], '%s+' %(c))
        keep_flag = raw_input('Keep selection? (y), n ') 
        if keep_flag == 'n':
            ln[0].set_visible(False)
            pyplot.draw()

    return coords[0][0]
    
def select_extraction_box_size(fig1, extrlocy, c, region_to_mark = 'spectrum'): #you should never enter this function if mode = passive
    '''Interactively select the box size of an extraction or background region '''
    print 'Select extraction box for %s' %(region_to_mark)
    keep_flag = 'n'
    while keep_flag == 'n':
        coords = fig1.ginput(n = 2, timeout = -1)
        half_box = min(abs(coords[0][0] - extrlocy), abs(coords[1][0] - extrlocy))
        lb = pyplot.plot([extrlocy - half_box, extrlocy + half_box],[coords[0][1], coords[0][1]] , '%s|-' %(c), markersize = 10)
        pyplot.draw()
        keep_flag = raw_input('Keep points? (y), n ')
        if keep_flag == 'n':
            lb[0].set_visible(False)
            pyplot.draw()
    return half_box*2

def extract_spectrum(input, extrlocy, extract_box_size, background_loc1, background_size1, background_loc2, background_size2, c, backcorr_option, bksmode_option):
    '''Extract the spectrum using the calstis task x1d '''
    if os.path.exists(os.path.join(os.getcwd(), input.filename.replace('.fits', '_loc%i.fits' %(int(extrlocy))))):
        os.remove(os.path.join(os.getcwd(), input.filename.replace('.fits', '_loc%i.fits' %(int(extrlocy)))))
        
    #iraf.stsdas.hst_calib.stis.x1d(input.filename.replace('.fits', 'sub.fits'), output = input.filename.replace('.fits', '_loc%i.fits' %(int(extrlocy))), \
    #                                  a2center = extrlocy + 1, extrsize = extract_box_size, maxsrch = 0, bk1offst = background_loc1 - extrlocy, \
    #                                   bk2offst = background_loc2 - extrlocy, bk1size = background_size1, bk2size = background_size2, backcorr = backcorr_option, bksmode = bksmode_option)  #a2center is 1 indexed
    x1d(input.filename.replace('.fits', 'sub.fits'), output = input.filename.replace('.fits', '_loc%i.fits' %(int(extrlocy))), \
                                      a2center = extrlocy + 1, extrsize = extract_box_size, maxsrch = 0, bk1offst = background_loc1 - extrlocy, \
                                       bk2offst = background_loc2 - extrlocy, bk1size = background_size1, bk2size = background_size2, backcorr = backcorr_option, bksmode = bksmode_option)
    add_background_into_x1d(input.filename, fit, extrlocy)

def add_background_into_x1d(filename, fit, extrlocy):
    '''Add the fitted background back into the gross and background columns of the extracted spectrum'''
    ofile = pyfits.open(os.path.join(os.getcwd(), filename.replace('.fits', '_loc%i.fits' %(int(extrlocy)))), mode = 'update')
    exptime = ofile[0].header['texptime']
    ofile[1].data['gross'][:] = ofile[1].data['gross'][:] + fit[int(round(ofile[1].data['a2center'])) - 1]/exptime #a2center is 1 indexed, units are cts/s in the x1d, cts in the flt
    ofile[1].data['background'][:] = ofile[1].data['background'][:] + fit[int(round(ofile[1].data['a2center'])) - 1]/exptime #a2center is 1 indexed
    ofile.flush()
    ofile.close()

def confirm_extraction_location(filename, img, extrlocy):
    fig3 = pyplot.figure(3)
    ax3 = fig3.add_subplot(1, 1, 1)
    ax3.imshow(np.log10(img), interpolation = 'nearest', origin = 'lower', cmap = 'gray')
    tbdata =  pyfits.getdata(os.path.join(os.getcwd(), filename.replace('.fits', '_loc%i.fits' %(int(extrlocy)))),  1)
    extrlocy = tbdata['extrlocy'].ravel() - 1.0
    extrsize = tbdata['extrsize'].ravel() 
    bk1size = tbdata['bk1size'].ravel()
    bk2size = tbdata['bk2size'].ravel()
    bk1offset = tbdata['bk1offst'].ravel()
    bk2offset = tbdata['bk2offst'].ravel()
    pix = np.arange(len(extrlocy))
    pyplot.plot(pix, extrlocy, 'r', lw = 3)
    pyplot.plot(pix, extrlocy - 0.5*extrsize, 'r--', lw = 2)
    pyplot.plot(pix, extrlocy + 0.5*extrsize, 'r--', lw = 2)
    pyplot.plot(pix, extrlocy + bk1offset, color = '#0022FF', lw = 2)
    pyplot.plot(pix, extrlocy + bk1offset + 0.5*bk1size, color = '#0022FF', ls = '--', lw = 2)
    pyplot.plot(pix, extrlocy + bk1offset - 0.5 * bk1size, color = '#0022FF', ls = '--', lw = 2)
    pyplot.plot(pix, extrlocy + bk2offset, color = '#0022FF', lw = 2)
    pyplot.plot(pix, extrlocy + bk2offset + 0.5*bk2size, color = '#0022FF', ls = '--', lw = 2)
    pyplot.plot(pix, extrlocy + bk2offset - 0.5*bk2size, color = '#0022FF', ls = '--', lw = 2)
    raw_input('Press enter to close and continue')
    pyplot.close(fig3)
    
def plot_final_spectrum(filename, extrlocy):
    fig4 = pyplot.figure(4)
    ax4 = fig4.add_subplot(1, 1, 1)
    tbdata = pyfits.getdata(os.path.join(os.getcwd(), filename.replace('.fits', '_loc%i.fits' %(int(extrlocy)))),  1)
    ax4.plot(tbdata['wavelength'][0], tbdata['flux'][0])
    ax4.set_xlabel('Wavelength')
    ax4.set_ylabel('Flux')
    ax4.set_title('Spectrum for %s @ %i' %(filename, int(extrlocy)))
    raw_input('Press enter to close and continue')
    pyplot.close(fig4)

def interact_w_user(mode, default = None, message = None):
    if (mode == 'interactive'): 
        if message is None:
            print 'You must enter a message if mode is interactive'
            sys.exit() 
        else:
            result = raw_input(message)
            return result
    else:
        return default

class input_object:
    def __init__(self, filename, mode):
        self.filename = filename
        self.mode = mode

if __name__ == "__main__":

    #Define colors for extracting more than one spectrum
    colors = ['r', 'g', 'c', 'k', 'm']

    pyplot.ion()  #turn plotting on 
    if os.path.exists('/grp/hst/cdbs/oref'):
        os.environ['oref'] = '/grp/hst/cdbs/oref/' #set oref environment variable to point to reference file location
        os.environ['myref'] = '/user/bostroem/science/cte/2012_04/reffiles/'
    else:
        os.environ['oref'] = '/Users/bostroem/science/oref/'
        os.environ['myref'] = '/Users/bostroem/science/cte/2012_04/reffiles/'

    parser = OptionParser()
    parser.add_option('--backcorr', dest = 'backcorr', help = 'Enter perform (default) or omit to perform or omit the background subtraction in CalSTIS x1d', default = 'perform')
    parser.add_option('--ncol', dest = 'num_cols', type = 'float', help = 'Number of columns summed when examining the cross-dispersion profile', default = 50)
    parser.add_option('--backsmooth', dest = 'bksmode', help = 'Background smoothing mode: off, median, average', default = 'off')
    parser.add_option('--mode', dest = 'mode', help = 'Run program interactively or not - options: interactive, passive', default = 'interactive')
    (options, args) = parser.parse_args()
    
    print options.bksmode
    #read in command line arguments
    filename = sys.argv[1]
    try:
        ext = int(sys.argv[2])
    except:
        ext = 1

    input = input_object(filename, options.mode)

    cenwave = pyfits.getval(filename, 'cenwave', 0)
    if os.path.exists(os.path.join(os.getcwd(), filename.replace('.fits','_c%i_background_regions.txt'%(cenwave)))):
        find_background_flag = interact_w_user(input.mode, default = 'n', message = 'Select the background regions? y, (n) ')
    else:
        find_background_flag = 'y'

    #set up figure and plot cross-dispersion profile
    fig1 = pyplot.figure(1)
    ax1 = fig1.add_subplot(1, 1, 1)
    ax1.set_xlabel('Pixel Number')
    ax1.set_ylabel('Total Counts Summed in the Dispersion Direction')
    #pdb.set_trace()
    img = pyfits.getdata(filename, ext)
    
    collapsed_img = collapse_spectrum(img, options.num_cols)
    pix_num = np.arange(len(collapsed_img))
    ax1.plot(pix_num, collapsed_img)
    ax1.set_ylim(np.min(collapsed_img), np.max(collapsed_img))
    if find_background_flag == 'y':  #interactvely define background
        #Display 2D image
        disp_2d = interact_w_user(input.mode, default = 'n', message = 'Would you like to display the 2D image? y, (n) ')
        if disp_2d == 'y':
            fig2 = pyplot.figure(2)
            ax2 = fig2.add_subplot(1,1,1)
            ax2.imshow(np.rot90(np.log10(img)), interpolation = 'nearest', origin = 'lower')
        interact_w_user(input.mode, message = 'Zoom in on desired spectrum location. Press Enter when finished')
        
        #Define background regions in cross dispersion profile
        fig1, ax1, background_pix, background_collapsed = get_background_locations(input, fig1, ax1, pix_num, collapsed_img, cenwave)

    else:  #read background regions from file
        coords = np.genfromtxt(filename.replace('.fits','_c%i_background_regions.txt'%(cenwave)))
        background_collapsed = np.empty((0,))
        background_pix = np.empty((0,))
        for region in coords:
            background_pix = np.append(background_pix, pix_num[int(region[0]):int(region[1])])
            background_collapsed = np.append(background_collapsed, collapsed_img[int(region[0]):int(region[1])])
        sort_indx = np.argsort(background_pix)
        background_pix = background_pix[sort_indx]
        background_collapsed = background_collapsed[sort_indx]

    #Get limits now to scale plot when plotting the background subtracted data (things can go crazy at the edges)
    xlims = ax1.get_xlim()
    ylims = ax1.get_ylim()
    
    #Fit a polynomial to the background
    ax1_temp, fit = fit_background(input, fig1, ax1, background_pix, background_collapsed, pix_num)
    while ax1_temp is None:
        fig1, ax1, background_pix, background_collapsed = get_background_locations(input, fig1, ax1, pix_num, collapsed_img, cenwave, first_time = False)
        ax1_temp, fit = fit_background(input, fig1, ax1, background_pix, background_collapsed, pix_num)
    ax1 = ax1_temp
    #Subtract the background (in the cross dispersion direction) from the image
    subtracted_img = subtract_2D_cross_disp_background(img, fit)
    #Write a temporary file of the subtracted image
    write_temp_subtracted_file(input, subtracted_img)
    #collapse the subtracted image in the cross-dispersion direction
    sub_collapsed_img = collapse_spectrum(subtracted_img, options.num_cols)
    ax1.cla()
    ax1.plot(pix_num, sub_collapsed_img) #Plot the cross dispersion profile of the subtracted image
    ax1.axhline(0, ls = '--', color = 'k')
    another_spectrum_flag = 'y'  #this flag denotes wanting to extract more than one spectrum from a subtracted image
    i = 0
    while another_spectrum_flag != 'n':  #for each spectrum
        c = colors[i%len(colors)] #Set color for spectrum extraction regions
        ax1.set_xlim(xlims[0], xlims[1])
        ax1.set_ylim(ylims[0], ylims[1])
        pyplot.draw()
        #print '---------Identify spectrum location---------'
        extrlocy = select_extraction_location(fig1, c)
        #print '---------Identify extraction box---------'
        extract_box_size = select_extraction_box_size(fig1, extrlocy, c)
        #print '---------Identify left background center---------'
        background_loc1 = select_extraction_location(fig1, c, region_to_mark = 'center of first background region')
        #print '---------Idnetify left background box---------'
        background_size1 = select_extraction_box_size(fig1, background_loc1, c, region_to_mark = 'first background region')
        #print '---------Identify right background center---------'
        background_loc2 = select_extraction_location(fig1, c, region_to_mark = 'second background region')
        #print '---------Identify right background box---------'
        background_size2 = select_extraction_box_size(fig1, background_loc2, c, region_to_mark = 'second background region')

        extract_spectrum(input, extrlocy, extract_box_size, background_loc1, background_size1, background_loc2, background_size2, c, options.backcorr, options.bksmode)
        confirm_flag = interact_w_user(input.mode, message = 'Would you like to confirm the location of your extraction on a 2D image? (y), n ', default = 'n')
        if confirm_flag != 'n':
            confirm_extraction_location(input.filename, img, extrlocy)
        final_spec_flag = interact_w_user(input.mode, message = 'Would you like to plot the final spectrum? (y), n ', default = 'n')
        if final_spec_flag != 'n':
            plot_final_spectrum(filename, extrlocy)
        another_spectrum_flag = interact_w_user(input.mode, message = 'Extract another spectrum? ', default = 'n')
        i = i + 1
    os.remove(input.filename.replace('.fits', 'sub.fits'))



#Add log file which records degree polynomial fit, extraction locations and box sizes
#To make mode = passive work, have to add way to select extraction and background location and regions
