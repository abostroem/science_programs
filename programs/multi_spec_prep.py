import pylab
import pyfits
import numpy as np
import os
import glob
import pdb

def show_image(filename, disp_log = True, ext = 0):
    img = pyfits.getdata(filename, ext)
    fig = pylab.figure('1', figsize = [25, 20])
    ax = fig.add_subplot(1,1,1)
    if disp_log:
        ax.imshow(np.log10(img), interpolation = 'nearest', origin = 'lower')
    else:
        ax.imshow(img, interpolation = 'nearest', origin = 'lower')
    ax.axvline(np.shape(img)[1]/2, color = 'grey', ls = '--', lw = 3)
    spec_loc = get_spectrum_locations(fig)
    for xy_coor in spec_loc:
        pylab.plot(xy_coor[0], xy_coor[1], '.', color = 'r', markersize = 8, markeredgecolor = 'k')
    raw_input('Press enter to continue')
    pylab.close()
    return spec_loc


def get_spectrum_locations(fig):
    print 'Left click (or press any key) to select a point'
    print 'Right click (or press backspace or delete) to pop the last point input'
    print 'Middle click (or press enter) to terminate input'

    #timeout = -1 --> will never timeout
    #n = -1 --> no set number of clicks
    spec_loc = np.array(fig.ginput(-1, timeout = -1))
    return spec_loc

def write_output(out_filename, spec_loc, idir = '/user/bostroem/science/multi_spec_files/'):
    ofile = open(os.path.join(idir, out_filename), 'w')
    ofile.write('id\t x_center\ty_center\tRA\tDEC\tSED\te\n')
    for i, xy_coor in enumerate(spec_loc): 
        #E(B-V) recommended by Jesus Maiz Apellaniz
        ofile.write('%i\t\t%f\t%f\t%i\t%i\t\t0.5\n' %(i, xy_coor[0], xy_coor[1], 0, 0))
    ofile.close()




if __name__ == "__main__":
    idir = '../12465_otfr20120425/mama/'
    os.chdir(idir)
    flist = glob.glob('R136*combined_img.fits')
    for ifile in flist:
        spec_loc = show_image(ifile)
        write_output('%s_phot.dat' %(ifile[0:8]), spec_loc)
                 
              


