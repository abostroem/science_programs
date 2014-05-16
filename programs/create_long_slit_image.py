import pyfits
import pylab
import numpy as np
import glob
import os
import pdb
import scipy
import sys
import math
from scipy import ndimage
from scipy.ndimage import interpolation
from scipy.ndimage.interpolation import rotate
import matplotlib
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap

def collapse_longslit(spec, extract_loc_list = None, plate_scale = 0.0246, slit_width = 0.2):
    '''
    ########################################################################################################################
    #This function takes a 2D image (from an FLT or X2D file and collapses it in the dispersion direction)
    #It expands the 1D array to be the width of the slit and returns the new 2D 'slit' image
    #input:
    #    spec: a 2D image which will be collapsed
    #    plate_scale: the plate scale in arcsec/pix. Default to STIS G140L plate scale
    #    slit_width: the width of the slit used in arcsec
    #Output:
    #    spec_tile: an 'image' of the long slit
    #Called from:
    #    create_all_tile_files
    ########################################################################################################################
    '''
    one_d_array = np.sum(spec, axis = 1)
    if extract_loc_list:
        print np.max(one_d_array)
        for yloc in extract_loc_list:
            try:
                one_d_array[int(math.floor(yloc))] = 0.0
                one_d_array[int(math.ceil(yloc))] = 0.0
            except:
                pass
    #0.0246 arcsec/pix, 0.2 arcsec slit = 8.1 pixels
    slit_width_pix = round(slit_width/plate_scale)
    spec_tile = np.transpose(np.tile(one_d_array, (slit_width_pix,1)))
    return spec_tile, round(slit_width/plate_scale)

def write_fits_tile(spec_tile, filename):
    '''
    ########################################################################################################################
    #This program writes the spec_tile output from collapse_longslit to a fits file
    #Input:
    #    spec_tile: 2D image to write
    #    filename: name of file to write to (without extension)
    #Output:
    #    writes a fits file
    #Called from:
    #    create_all_tile_files
    #    display_rotated_image
    ########################################################################################################################
    '''
    print filename
    
    hdu = pyfits.PrimaryHDU(spec_tile)
    hdulist = pyfits.HDUList([hdu])
    #if os.path.exists('/user/bostroem/science/12465_otfr20120425/mama/%s.fits' %(filename)):
    #    os.remove('/user/bostroem/science/12465_otfr20120425/mama/%s.fits' %(filename))
    hdulist.writeto(filename+'.fits', clobber = True)

def create_all_tile_files(path_name, flist, extract_filename = None, overwrite = True):
    '''
    ########################################################################################################################
    #This is the main function to call to write all fits files of the slit images
    #Input:
    #    path_name: path where original data (and output data) live
    #    overwrite: should the fits files be overwritten? default is to overwrite
    #Ouptut:
    #    writes fits files for all of the slit images
    #    returns the floor of the slit width in pixels
    #Calls to:
    #    collapse_longslit
    #     write_fits_tile
    #Called from:
    #    __main__
    ########################################################################################################################
    '''
    os.chdir(path_name)
    opt_elem = pyfits.getval(flist[0], 'opt_elem', 0)
    cenwave = pyfits.getval(flist[0], 'cenwave', 0)
    if extract_filename:
        extract_dict = make_extraction_dictionary(extract_filename)
    for ifile in flist:
        print ifile
        spec = pyfits.getdata(ifile, 0)
        if extract_filename:
            extract_loc_list = extract_dict[ifile[0:8]]
        if opt_elem == 'G140L':
            if extract_filename:
                spec_tile, slit_size = collapse_longslit(spec, extract_loc_list = extract_loc_list, plate_scale = 0.0246, slit_width = 0.2)
            else:
                spec_tile, slit_size = collapse_longslit(spec, plate_scale = 0.0246, slit_width = 0.2)
        elif opt_elem == 'G430M':
            spec_tile, slit_size = collapse_longslit(spec, plate_scale = 0.05078, slit_width = 0.2)
        else:
            print 'Grating %s is not recognized' %(opt_elem)
            sys.exit()
        if not os.path.exists(path_name+'/'+ifile[0:8]+'_collapse_flt.fits'):
            write_fits_tile(spec_tile, ifile[0:8]+'_collapse_flt')
        else:
            if overwrite == True:
                #os.remove(path_name+'/'+ ifile[0:8]+'_collapse_flt.fits')
                write_fits_tile(spec_tile, ifile[0:8]+'_%i_collapse_flt' %(cenwave))
            else:
                print 'file already exists, not overwriting'
    return math.floor(slit_size)


def combine_images(flist, slit_size):
    '''
    ########################################################################################################################
    #This function takes a list of slit image files and the slit width and combines them into an image of multiple slits on the sky
    #Inputs:
    #    flist: list of slit image files
    #    slit_size: the width of the slit in pixels
    #Output:
    #    combined_image: an image of all of the slit images combined in a row
    #Called from:
    #    create_image
    ########################################################################################################################
    '''
    slit_size = int(slit_size)
    img = pyfits.getdata(flist[0], 0)
    combined_image = np.empty((np.shape(img)[0], slit_size*len(flist)))
    #combined_image = np.empty((1026, (slit_size+1.0)*len(flist)+1.0))
    for i, ifile in enumerate(flist):
        img = pyfits.getdata(ifile, 0)

        combined_image[:,slit_size*i:slit_size*(i+1)] = img
        #combined_image[1:1026,((slit_size+1)*i+1:(slit_size+1)*(i+1)] = img
    hdu = pyfits.PrimaryHDU(combined_image)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto('long_slit_image.fits', clobber = True)
    return combined_image


def mark_pix_loc(xdim, ydim, angle_rad, ax1):
    #Mark the lower edge
    x_orig = np.arange(ydim)
    x_orig_rev = x_orig[::-1]
    x_actual = -math.cos(angle_rad)*x_orig + np.max(math.cos(angle_rad)*x_orig)
    y_actual = math.sin(angle_rad)*x_orig
    x_actual_rev = x_actual[::-1]
    y_actual_rev = y_actual[::-1]
    ax1.plot(x_actual_rev[::50], y_actual_rev[::50], 'k|', markersize = 3, markeredgewidth = 2)
    for xo, xa, ya in zip(x_orig[::50], x_actual_rev[::50], y_actual_rev[::50]):
        ax1.text(xa, ya - 20, str(xo), fontsize = 8, ha = 'center')
    #Mark the upper edge
    x_actual_u = x_actual + xdim*math.sin(angle_rad)
    y_actual_u = y_actual + xdim * math.cos(angle_rad)
    x_actual_u_rev = x_actual_u[::-1]
    y_actual_u_rev = y_actual_u[::-1]
    ax1.plot(x_actual_u_rev[::50], y_actual_u_rev[::50], 'k|', markersize = 3, markeredgewidth = 2)
    for xo, xa, ya in zip(x_orig[::50], x_actual_u_rev[::50], y_actual_u_rev[::50]):
        ax1.text(xa, ya + 10, str(xo), fontsize = 8, ha = 'center')    
    return ax1

def mark_boundaries(flist, slit_size, combined_img, rot_img, ax1, img_angle, ref_angle, fudge_factor, target_font_size):
    flist_r = flist[::-1]

    angle =  (ref_angle - img_angle) - 90
    angle_rad = angle/180.0 * math.pi
    ydim, xdim = np.shape(combined_img)
    x = np.arange( (ydim) * math.cos(angle_rad))
    
    slope = -math.tan(angle_rad)
    for indx, i in enumerate(np.arange(xdim/int(slit_size)+1)):                             
        intercept = (ydim) * math.sin(angle_rad) + ((i)*slit_size/math.cos(angle_rad))
        xnew = x+((i)*math.sin(angle_rad)*slit_size)
        y = slope*xnew+intercept
        ax1.plot(xnew, y, 'k')
        try:
            if np.mod(i, 2) == 1:
                ax1.text(math.ceil(xnew[-1]+1), math.ceil(y[-1]), '%s' %(pyfits.getval(flist_r[indx].replace('collapse_flt', 'combined_img'), 'targname', 0)) , fontsize = target_font_size, ha = 'left', va = 'bottom')
                ax1.plot(np.arange(50) +  math.ceil(xnew[-1]), np.ones(50) * math.ceil(y[-1]) - 1, 'k')
                ax1.plot((-np.arange(50)) +  math.ceil(xnew[0]), np.ones(50) * math.ceil(y[0]) - 1, 'k')
            else:
                ax1.text(math.ceil(xnew[0]-1), math.ceil(y[0]), '%s' %(pyfits.getval(flist_r[indx].replace('collapse_flt', 'combined_img'), 'targname', 0)) , fontsize = target_font_size, ha = 'right', va = 'bottom')
                ax1.plot((-np.arange(50)) +  math.ceil(xnew[0]) - 1, np.ones(50) * math.ceil(y[0]) - 1, 'k')
                ax1.plot(np.arange(50) +  math.ceil(xnew[-1]) - 1, np.ones(50) * math.ceil(y[-1]) - 1, 'k')
        except:
            pass
    ax1 = mark_pix_loc(xdim, ydim, angle_rad, ax1)
    #pdb.set_trace()
    return ax1

def make_custom_colormap():
    cbar = pylab.get_cmap()
    cdict = cbar._segmentdata
    for key in cdict.keys():
        list1 = list(cdict[key])
        save_val = list(list1.pop(0))
        save_val[0] = 0.01

        list1.insert(0, tuple(save_val))
        list1.insert(0, (0.0, 1.0, 1.0))

        cdict[key] = tuple(list1)
    new_cmap = matplotlib.colors.LinearSegmentedColormap('colormap', cdict)
    #pdb.set_trace()
    return new_cmap
            
        
        
    

def rotate_image(img, img_angle, ref_angle, fudge_factor = 0):
    '''
    ########################################################################################################################
    #This function rotates the multi-slit image
    #Inputs:
    #    img: the multi-slit image array output by combine_images
    #    img_angle: the ORIENTAT angle from the slit image data
    #    ref_angle: the ORIENTAT angle from the WFC3 image that you are trying to match
    #    fudge_factor: An additional angle determined manually by the user to make the multi-slit image match the WFC3 image
    #Output:
    #    rot_img: the rotated image array
    #Called from:
    #    display_rotated_image
    #    display_rotated_image_and_wfc3_image
    ########################################################################################################################
    '''
    #angle =ref_angle - img_angle + fudge_factor
    angle = 360.0 - (ref_angle - img_angle)
    print 'rotating image %f deg E of N' %(angle)
    
    pdb.set_trace()
    rot_img = rotate(img, -angle) #rotate rotated CW (W from N), orientat is measured CCW (E from N)
    return rot_img

def display_rotated_image_and_wfc3_image(combined_image1, flist, wfc3_image, target_font_size, ff = -8.1, log = True, save_filename = 'rotated_img', cmap1 = 'jet', clim1 = None, clim2 = None, save = False, ax1_title = None, ax2_title = None):
    '''
    ########################################################################################################################
    #This function displays the rotated multi-slit image next to the WFC3 image
    #Inputs:
    #    combined_image1: the multi-slit image array
    #    wfc3_image: the file name of the WFC3 image
    #    ff: fudge factor used for additional rotation in rotate_image; default = -8.1
    #    log: display the log of the image; default = True
    #    save_filename: save images to this filename; default = 'rotated_img'
    #    cmap1: Color map to use; default = None - use default matplotlib colorbar
    #    clim1: lower contrast limit; default = None - use default matplotlib clim
    #    clim2: upper contrast limit; default = None - use default matplotlib clim
    #    save: switch to enable user to save the file (to save_filename); default = False
    #Output:
    #    if the save keyword is set then the images will be displayed
    #Calls to:
    #    rotate_image
    #Called from:
    #    create_image
    ########################################################################################################################
    '''
    wfc3_img  = pyfits.getdata(wfc3_image, 0)
    rot_img = rotate_image(combined_image1, 64.0072, 166.002207094, fudge_factor = ff)
    fig = pylab.figure(figsize = [30, 20])
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    #new_colormap = make_custom_colormap()
    if not cmap1: cmap1 = 'jet'
    #pdb.set_trace()
    new_colormap = getattr(matplotlib.cm, cmap1)
    new_colormap1 = getattr(matplotlib.cm, 'jet')
    norm1 = colors.Normalize(vmin = np.min(np.log10(wfc3_img)[np.isfinite(np.log10(wfc3_img))]) + 0.01, vmax = np.max(np.log10(wfc3_img)[np.isfinite(np.log10(wfc3_img))]))
    new_colormap.set_under('white')
    #pdb.set_trace()
    ax2.imshow(np.log10(wfc3_img), interpolation = 'nearest', cmap = new_colormap1, norm = norm1)

    if log:
        rot_img_log = np.log10(rot_img)
        nan_indx = np.isnan(rot_img_log)
        inf_indx = np.isinf(rot_img_log)
        neg_inf_indx = np.isneginf(rot_img_log)
        rot_img_log[nan_indx] = 0
        rot_img_log[inf_indx] = 0
        rot_img_log[neg_inf_indx] = 0
        rot_img_log[nan_indx] = np.min(rot_img_log) - 1
        rot_img_log[inf_indx] = np.min(rot_img_log) - 1
        rot_img_log[neg_inf_indx] = np.min(rot_img_log) - 1
        #norm2 = colors.Normalize(vmin = np.min(rot_img_log[np.isfinite(rot_img_log)]) + 0.01, vmax = np.max(rot_img_log[np.isfinite(rot_img_log)]))
        norm2 = colors.Normalize(vmin = 0, vmax = np.max(rot_img_log[np.isfinite(rot_img_log)]))
        cax = ax1.imshow(rot_img_log, interpolation = 'nearest', cmap = new_colormap, norm = norm2)
    else:
        norm2 = colors.Normalize(vmin = np.min(rot_img) + 0.01, vmax = np.max(rot_img))
        cax = ax1.imshow(rot_img, interpolation = 'nearest')
    ax2.set_xlim(1300, 2000)
    ax2.set_ylim(1500, 2100)
    ax1.set_xlim(-100, 1100)
    ax1.set_ylim(-20, 490)
    fig.colorbar(cax)
    #if cmap1:
    #    cax.set_cmap(cmap1)
    #pdb.set_trace()
    if clim1:
        cax.set_clim(clim1, clim2)
    ax1 = mark_boundaries(flist, slit_size, combined_image1, rot_img, ax1, 64.0072, 166.002207094, fudge_factor = ff, target_font_size = target_font_size)
    if ax1_title:
        ax1.set_title(ax1_title)
    if ax2_title:
        ax2.set_title(ax2_title)
    pdb.set_trace()
    if save:
        pylab.savefig(save_filename+'.pdf')
    

def create_image(path_name, flist, slit_size,  start_indx = 0, skip_indx = 1, log = True, save_filename = 'rotated_image', cmap1 = None, clim1 = None, clim2 = None,save = False, ax1_title = None, ax2_title = None,  target_font_size = 6):
    '''
    ########################################################################################################################
    #This function takes a set of slit image fits files and creates a single image out of them and displays them with a WFC3
    #image
    #Inputs:
    #    path_name: directory where slit image fits files exist
    #    glob_str1: universal string to select slit image files using glob
    #    slit_size: size of the slit in pixels
    #    glob_str2: second optional  universal string to select slit image files using glob; default = None
    #    start_indx: start using files at this index; default = 0
    #    skip_indx: increment file list by this number. This is so you can select every other file (for instance); default = 1
    #    log: display the log of the image; default = True
    #    save_filename: save images to this filename; default = 'rotated_img'
    #    cmap1: Color map to use; default = None - use default matplotlib colorbar
    #    clim1: lower contrast limit; default = None - use default matplotlib clim
    #    clim2: upper contrast limit; default = None - use default matplotlib clim
    #    save: switch to enable user to save the file (to save_filename); default = False
    #Output:
    #    combined_image1: multi-slit image array
    #Calls to:
    #    combine_images
    #    display_rotated_image_and_wfc3_image
    #Called by:  __main__
    ########################################################################################################################
    '''
    os.chdir(path_name)
    combined_image1 = combine_images(flist, slit_size)
    #display_rotated_image(combined_image1)
    display_rotated_image_and_wfc3_image(combined_image1, flist[start_indx::skip_indx], '/user/bostroem/science/images/f336_med_drz_sci.fits', ff = 150, log = log, save_filename = save_filename, cmap1 = cmap1, clim1 = clim1, clim2 = clim2, save = save, ax1_title = ax1_title, ax2_title = ax2_title,  target_font_size = target_font_size)
    pdb.set_trace()
    return combined_image1


###################TEST FUNCTIONS####################

def display_joined_images(combined_image1, combined_image2): 
    '''
    ########################################################################################################################
    #The data that we have are dithered along the slit. This creates and displays an image from each dither locations (so you get
    #essentially 2 identical images
    #Inputs:
    #    combined_image1: first image array to display
    #    combined_image2: second image array to display
    ########################################################################################################################  
    '''
    fig = pylab.figure()
    ax1 = fig.add_subplot(1,2,1)
    ax1.imshow(np.log10(combined_image1), interpolation = 'nearest')
    ax2 = fig.add_subplot(1,2,2)
    ax2.imshow(np.log10(combined_image2), interpolation = 'nearest')
    pdb.set_trace()

def display_rotated_image(combined_image1, ff = -8.1):
    '''
    ########################################################################################################################
    #Display the muli-slit image rotated by the difference of the WFC3 image and the slit images and a fudge factor
    #Inputs:
    #    combined_image1: the multi-slit array to combined
    #    ff: the fudge factor used for additional rotation in rotate_image
    #Calls to:
    #    rotate_image
    #    write_fits_tile
    #Called from:
    ########################################################################################################################
    '''
    rot_img = rotate_image(combined_image1, 64.0072, 166.002207094, fudge_factor = ff)
    fig = pylab.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.imshow(np.nan_to_num(np.log10(rot_img)), interpolation = 'nearest')
    write_fits_tile(rot_img, 'rotated_img')
    pdb.set_trace()

def make_extraction_dictionary(extraction_filename):
    extract_dict = {}
    data_struc = np.genfromtxt(extraction_filename, names = ['target', 'extrlocy'], skiprows = 0, usecols = [0, 1], dtype = "S8, f8")
    for targ in set(data_struc['target']):
        extract_dict[targ] = []
    for targ, loc in zip(data_struc['target'], data_struc['extrlocy']):
        extract_dict[targ].append(loc)
    return extract_dict

if __name__ == "__main__":

    #set plotting parameters
    pylab.rcParams['figure.subplot.left'] = 0.035
    pylab.rcParams['figure.subplot.right'] =  0.99
    pylab.rcParams['figure.subplot.wspace'] = 0.07
    flist = ['R136-SE9_combined_img.fits', 'R136-SE8_combined_img.fits', 'R136-SE7_combined_img.fits', 'R136-SE6_combined_img.fits', 'R136-SE5_combined_img.fits', 'R136-SE4_combined_img.fits', 'R136-SE3_combined_img.fits', 'R136-SE2_combined_img.fits', 'R136-SE1_combined_img.fits', 'R136-NW1_combined_img.fits', 'R136-NW2_combined_img.fits', 'R136-NW3_combined_img.fits', 'R136-NW4_combined_img.fits', 'R136-NW5_combined_img.fits', 'R136-NW6_combined_img.fits', 'R136-NW7_combined_img.fits', 'R136-NW8_combined_img.fits']
    flist2 = []
    for ifile in flist:
        flist2.append(ifile.replace('combined_img', 'collapse_flt'))

    ############################
    #Create MAMA Long Slit Image
    ############################
    #Create slit image fits files
    #slit_size = create_all_tile_files('/user/bostroem/science/12465_otfr20120425/mama', flist, extract_filename = 'extraction_locations_fuv.dat') #,  overwrite = False)   #for the FUV-MAMA
    #slit_size = create_all_tile_files('/user/bostroem/science/12465_otfr20120425/mama', flist) #,  overwrite = False)   #for the FUV-MAMA


    #Test dithered images
    #img1 = create_image('/user/bostroem/science/12465_otfr20120425/mama', 'obrc04???_collapse_flt.fits', glob_str2 = 'obrc05???_collapse_flt.fits', skip_indx = 2)
    #img2 = create_image('/user/bostroem/science/12465_otfr20120425/mama', 'obrc04???_collapse_flt.fits', glob_str2 = 'obrc05???_collapse_flt.fits', start_indx = 1, skip_indx = 2)
    #display_joined_images(img1, img2)

    #Use slit image fits files to display multi-slit image and WFC3 image

    #MAMA
    #create_image('/user/bostroem/science/12465_otfr20120425/mama', flist2, slit_size, skip_indx = 1, \
    #                 save_filename = 'rotated_img_mama_c1425', save = True, \
    #                 ax1_title = 'G140L Long-Slit Images from Visits 4 and 5', ax2_title = 'WFC3 F336 Image')#, clim1 = -2, clim2 = 6.1)



    ############################
    #Create CCD Long Slit Image
    ############################
    
    #CCD
    ##Create slit image fits files
    flist = ['R136-NW8_3936_combined_img.fits', 'R136-NW7_3936_combined_img.fits', 'R136-NW6_3936_combined_img.fits', 'R136-NW5_3936_combined_img.fits', 'R136-NW4_3936_combined_img.fits', 'R136-NW3_3936_combined_img.fits', 'R136-NW2_3936_combined_img.fits', 'R136-NW1_3936_combined_img.fits', 'R136-SE1_3936_combined_img.fits', 'R136-SE2_3936_combined_img.fits', 'R136-SE3_3936_combined_img.fits', 'R136-SE4_3936_combined_img.fits', 'R136-SE5_3936_combined_img.fits', 'R136-SE6_3936_combined_img.fits', 'R136-SE7_3936_combined_img.fits', 'R136-SE8_3936_combined_img.fits', 'R136-SE9_3936_combined_img.fits']

    flist2 = []
    for ifile in flist:
        flist2.append(ifile.replace('combined_img', 'collapse_flt'))
    flist2.reverse()
    slit_size = create_all_tile_files('/user/bostroem/science/12465_otfr20121109/ccd', flist, extract_filename = 'extraction_locations_ccd.dat')    #For the CCD
    create_image('/user/bostroem/science/12465_otfr20121109/ccd', flist2, slit_size, skip_indx = 1,  \
                     ax1_title = 'CCD G430L c3936 Long-Slit Images from Visit 1',ax2_title = 'WFC3 F336 Image',  target_font_size = 2, \
                     save_filename = 'rotated_img_CCD_c3936', save = True)# , cmap1 = 'spectral', clim1 = -2, clim2 = 6.1)

    #create_image('/user/bostroem/science/12465_otfr20120425/ccd', flist2, slit_size, skip_indx = 1, \
    #                 save_filename = 'rotated_img_CCD_c4451', save = True, cmap1 = 'spectral', clim1 = -2, clim2 = 6.1, target_font_size = 2, \
    #                 ax1_title = 'CCD G430L c4451 Long-Slit Images from Visit 2', ax2_title = 'WFC3 F336 Image')

    #create_image('/user/bostroem/science/12465_otfr20121109/ccd', flist2, slit_size, skip_indx = 1, \
    #                 save_filename = 'rotated_img_CCD_c4706', save = True, cmap1 = 'spectral', clim1 = -2, clim2 = 6.1, target_font_size = 2, \
    #                 ax1_title = 'CCD G430L c4706 Long-Slit Images from Visit 3', ax2_title = 'WFC3 F336 Image')
