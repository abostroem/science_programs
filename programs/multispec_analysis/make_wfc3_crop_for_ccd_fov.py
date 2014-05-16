#---------------------------
#Make WFC3 Image with CCD field of view
#---------------------------
from astropy.io import fits
from matplotlib import pyplot
import numpy as np


def crop_wfc3_img(img_file, xmin, xmax, ymin, ymax, output_filename):
    img =  fits.getdata(img_file, 0)
    cropped_img = img[ymin:ymax, xmin:xmax]
    pyplot.figure()
    pyplot.imshow(cropped_img, interpolation = 'nearest', cmap = 'bone', vmin = 0, vmax = 200)
    pyplot.savefig('../../multispec/ccd_multispec/{}'.format(output_filename.replace('.fits', '.pdf')))
    hdu = fits.PrimaryHDU(cropped_img)
    hdu.header.set('Instrume', 'WFC3')
    hdu.header.set('Proposid', 11360)
    hdu.header.set('xlim', '({}, {})'.format(xmin, xmax))
    hdu.header.set('ylim', '({}, {})'.format(ymin, ymax))
    hdu.header.set('EXPTIME', fits.getval(img_file,'exptime', 0))
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(output_filename, clobber = True)

if __name__ == "__main__":
    #crop_wfc3_img('../../images/astrodrizzle_wfc3_ccd_platescale/63.5_rot/final_drz_sci.fits',  2376+17, 2444+17, 814+17, 1858+17, '../../multispec/ccd_multispec/f336w_63.5_cropped.fits')
    #crop_wfc3_img('../../images/astrodrizzle_wfc3_ccd_platescale/64.06_rot/final_drz_sci.fits' 2376, 2444,, 814, 1858, '../../multispec/ccd_multispec/f336w_64.06_cropped.fits')
    #crop_wfc3_img('../../images/astrodrizzle_wfc3_ccd_platescale/63.6_rot/final_drz_sci.fits', 2376 + 12, 2444 + 12, 814+14, 1858+14, '../../multispec/ccd_multispec/f336w_63.6_cropped.fits')
    #crop_wfc3_img('../../images/astrodrizzle_wfc3_ccd_platescale/63.80_rot/final_drz_sci.fits', 2376 + 8, 2444 + 8, 814+10, 1858+10, '../../multispec/ccd_multispec/f336w_63.8_cropped.fits')
    crop_wfc3_img('../../images/astrodrizzle_wfc3_ccd_platescale/63.70_rot/final_drz_sci.fits', 2376 + 11, 2444 + 11, 814+12, 1858+12, '../../multispec/ccd_multispec/f336w_63.7_cropped.fits')