import numpy as np
import pdb

def get_wfc3_cts(wfc3_x, wfc3_y, wfc3_mag):
    #read in the .mag files from DAOPHOT to get the SUM value for each star
    ofile = open(wfc3_mag, 'r')
    all_lines = ofile.readlines()
    cts = []
    for iline in all_lines:
        if iline[0] != '#':  #skip the comments
            if iline.split()[0] == 'f336w_crop.fits':   #identify the start of data for a new star
                row_num = 0
            else:
                row_num += 1
            if row_num == 1:  #Check that you are matching the correct star to the correct counts
                if (float(iline.split()[0]) != wfc3_x[len(cts)]) or (float(iline.split()[1]) != wfc3_y[len(cts)]):
                    print iline.split()[0:2]
                    pdb.set_trace()
            if row_num == 5: #append the counts
                cts.append(float(iline.split()[1]))
    return np.array(cts)


def get_wfc3_coords(wfc3_mag):
    ofile = open(wfc3_mag, 'r')
    all_lines = ofile.readlines()
    wfc3_x = []
    wfc3_y = []
    cts = []
    for iline in all_lines:
        if iline[0] != '#':  #skip the comments
            if iline.split()[0] == 'f336w_crop.fits[0]':   #identify the start of data for a new star
                row_num = 0
            else:
                row_num += 1
            if row_num == 1:  
                wfc3_x.append(float(iline.split()[0])) 
                wfc3_y.append(float(iline.split()[1]))
            if row_num == 4: #append the counts  #This number is affected by the number of apertures you look at
                cts.append(float(iline.split()[1]))
                #pdb.set_trace()
    return np.array(wfc3_x), np.array(wfc3_y), np.array(cts)


def get_slit_num(stis_x):
    #Get the slit number. The first slit is 1. Pixels should start at 1 at the center of the first pixel
    slit_num = np.int_(np.int_(stis_x - 1)/8.0 + 1)
    return slit_num


def calc_y_stis_wrt_img_center(wfc3_y, y_corr, slit_num):
    #Calculated the STIS location based on the WFC3 location, average offset, and drift correction
    #Calculate the drift correction
    delta_y = np.float_(slit_num) 
    visit04_indx = np.where(delta_y <= 12) #drift correction is visit dependent
    visit05_indx = np.where(delta_y > 12)
    delta_y[visit04_indx] = 0.09489561*slit_num[visit04_indx] - 0.40810853
    delta_y[visit05_indx] = 0.06298826*slit_num[visit05_indx] - 1.55922388

    stis_y_corr = wfc3_y + y_corr - delta_y  #dist centered on 0 --> mean (or median) is a pretty good estimate of the offset
    return stis_y_corr

def calc_x_offset_wrt_slit_center(wfc3_x, x_corr):
    #Calcaulate the STIS location based on the WFC3 location and average offset between the images
    stis_x_corr = wfc3_x + x_corr #dist centered on 0 --> mean (or median) is a pretty good estimate of the offset
    slit_num = get_slit_num(stis_x_corr)
    stis_x_corr_slit_centered = stis_x_corr - ((slit_num - 1.0)*8.0 + 4.5)
    #pdb.set_trace()
    return stis_x_corr_slit_centered, slit_num
    
def write_ds9_reg_file(x, y, text = [], filename = 'ds9.reg'):
    ofile = open(filename, 'w')
    ofile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \nimage \n')

    if len(text) != 0:
        for ix, iy, it in zip(x, y, text):
            ofile.write('point(%f,%f) # point=cross text={%s} \n' %(ix, iy, str(it)))
    else:  
        for ix, iy in zip(x, y):
            ofile.write('point(%f,%f) # point=cross \n' %(ix, iy))
    ofile.close()


def calc_wfc3_stis_shift():
    #set file names
    stis_mag = '/user/bostroem/science/images/astrodrizzle_wfc3/final/add_more_stars/phot_output_stis.mag'
    wfc3_mag = '/user/bostroem/science/images/astrodrizzle_wfc3/final/add_more_stars/phot_output_wfc3.mag'
    matched_tab = '/user/bostroem/science/images/astrodrizzle_wfc3/final/add_more_stars/pair_coords.txt'
    #read in data
    wfc3_x, wfc3_y, stis_x, stis_y = np.genfromtxt(matched_tab,unpack = True)
    x_corr =  np.median(stis_x - wfc3_x)
    y_corr = np.median(stis_y - wfc3_y)
    return x_corr, y_corr
    
    

def create_table_from_geofile_analysis():
    x_corr, y_corr = calc_wfc3_stis_shift()

    #perform calculations
    cts = get_wfc3_cts(wfc3_x, wfc3_y, wfc3_mag)
    stis_x_corr, slit_num = calc_x_offset_wrt_slit_center(wfc3_x, x_corr)
    stis_y_corr = calc_y_stis_wrt_img_center(wfc3_y, y_corr, slit_num)
    write_ds9_reg_file(wfc3_x + np.median(stis_x - wfc3_x), stis_y_corr, slit_num)
    write_output_file(wfc3_x, wfc3_y, cts, slit_num, stis_x_corr, stis_y_corr)

def filter_duplicates(wfc3_x, wfc3_y, wx, wy, i):
    r = np.sqrt((wfc3_x - wx)**2 + (wfc3_y - wy)**2)
    dup_indx = np.where(r[:i] < 1.0)
    if len(dup_indx[0]) == 0:  #No Duplicates
        return True
    else:
        return False


def write_output_file(wfc3_x, wfc3_y, cts, slit_num, stis_x_corr, stis_y_corr, filename = 'multispec_input.txt', maxy = 1036):
    #write output
    ofile = open(filename, 'w')
    ofile.write('%10s    %10s    %10s    %8s    %10s    %10s \n' %('wfc3_x', 'wfc3_y', 'wfc3_cts', 'slit_num', 'stis_x', 'stis_y'))  
    for i, wx, wy, wct, sn, sx, sy in zip(range(len(wfc3_x)), wfc3_x, wfc3_y, cts, slit_num, stis_x_corr, stis_y_corr):
        #print i
        #if (np.abs(sx) <= 4.0) & ((sy >= 0.0) & (sy <= maxy)) & ((sn >= 1) & (sn <= 17)):
        if ((sy >= 0.0) & (sy <= maxy)) & ((sn >= 1) & (sn <= 17)):
            if filter_duplicates(wfc3_x, wfc3_y, wx, wy, i) is True:
                ofile.write('%10.6f    %10.6f    %10.3f    %8i    %10.6f    %10.6f \n' %(wx, wy, wct, sn, sx, sy))
            else:
                print 'duplicate', wx, wy, sx, sy, sn
        else:
            #print 'out of bounds', wx, wy, sx, sy, sn
            if ((sn < 1) | (sn > 17)): print 'slit out of bound', wx, wy, sx, sy, sn
            if ((sy < 0.0) | (sy > maxy)): print 'y is out of bounds', wx, wy, sx, sy, sn
    ofile.close()

def create_ds9_reg_for_source_list(filename, x_corr):
    wfc3_x, wfc3_y, slit_num, stis_y = np.genfromtxt(filename, usecols = (0, 1, 3, 5), skiprows = 1, unpack = True)
    write_ds9_reg_file(wfc3_x, wfc3_y, filename = '/user/bostroem/science/images/astrodrizzle_wfc3/final/add_more_stars/wfc3_unfiltered_source_list2.reg')
    stis_x = wfc3_x + x_corr
    write_ds9_reg_file(wfc3_x, wfc3_y, filename = '/user/bostroem/science/images/astrodrizzle_wfc3/final/add_more_stars/wfc3_source_list.reg')
    write_ds9_reg_file(stis_x, stis_y, slit_num, filename = '/user/bostroem/science/images/astrodrizzle_wfc3/final/add_more_stars/stis_source_list.reg')

def create_table_from_source_list():
    wfc3_mag = '/user/bostroem/science/images/astrodrizzle_wfc3/final/add_more_stars/phot_output_wfc3_combined_sourcelist2.mag'
    wfc3_x, wfc3_y, cts = get_wfc3_coords(wfc3_mag)
    #wfc3_idl =  '/user/bostroem/science/images/astrodrizzle_wfc3/final/add_more_stars/idl_aper_output_test5.txt'
    #wfc3_x, wfc3_y, cts = np.genfromtxt(wfc3_idl, unpack = True)

    x_corr, y_corr = calc_wfc3_stis_shift()
    stis_x_corr, slit_num = calc_x_offset_wrt_slit_center(wfc3_x, x_corr)
    stis_y_corr = calc_y_stis_wrt_img_center(wfc3_y, y_corr, slit_num)
    write_output_file(wfc3_x, wfc3_y, cts, slit_num, stis_x_corr, stis_y_corr, filename = '/user/bostroem/science/images/astrodrizzle_wfc3/final/add_more_stars/multispec_input_source_list.txt')
    create_ds9_reg_for_source_list('multispec_input_source_list.txt', x_corr)


                                                    

if __name__ == "__main__":
    create_table_from_source_list()

