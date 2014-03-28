import glob
import pyfits
import numpy as np
from matplotlib import pyplot
import os
import pdb

def make_empty_array(flist):
    #(0.2 arcsec) / (0.05078 arcsec/pix)  = 3.938558487593541
    array_width = len(flist) * 4
    array_height = 1044
    long_slit_img = np.zeros((array_height, array_width))
    return long_slit_img

def make_slit_img(filename, width = 4):
    img = pyfits.getdata(filename, 1)
    xd_profile = np.median(img[:, 487:537], axis = 1)
    slit_img = xd_profile
    for i in range(width-1):
        slit_img = np.vstack((slit_img, xd_profile))
    return slit_img.T

if __name__ == "__main__":
    width = 4
    target_names = ['SE9', 'SE8', 'SE7', 'SE6', 'SE5', 'SE4', 'SE3', 'SE2', 'SE1', 'NW1', 'NW2', 'NW3', 'NW4', 'NW5', 'NW6', 'NW7', 'NW8']
    flist = [os.path.join('/Users/bostroem/science/12465_otfr20121109/ccd/', '%s_3936_combined_img.fits' %(targ)) for targ in target_names]

    long_slit_img = make_empty_array(flist)
    for i, ifile in enumerate(flist):
        slit_img = make_slit_img(ifile)
        orientat = pyfits.getval(ifile, 'orientat', 1)
        height = np.shape(slit_img)[0]
        if orientat > 0:  #right side up
            long_slit_img[0:height, i*width: (i+1)*width] = slit_img
        else: #upside down
            #pass
            long_slit_img[1044 - height:, i*width: (i+1)*width] = np.rot90(slit_img, k = 2)  #Figure this out more scientifically?
    im = pyplot.imshow(long_slit_img, interpolation = 'nearest', vmin = 0, vmax = 100, cmap = 'bone')

    #Array is 1044 while slit images are 1031 so this is a rotation by 90 and an offset of 13 pixels
    hdu = pyfits.PrimaryHDU(long_slit_img)
    hdu.header.set('Instrume', 'STIS')
    hdu.header.set('Proposid', 12465 )
    hdu.header.set('TEXPTIME', pyfits.getval(ifile, 'texptime', 0))
    hdulist = pyfits.HDUList([hdu])
    hdu.writeto('/Users/bostroem/science/multispec/ccd_multispec/long_slit_img_%i.fits' %(pyfits.getval(ifile, 'cenwave', 0)), clobber = True)