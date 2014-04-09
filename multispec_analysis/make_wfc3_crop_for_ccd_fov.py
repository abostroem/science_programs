#---------------------------
#Make WFC3 Image with CCD field of view
#---------------------------
from astropy.io import fits
from matplotlib import pyplot
import numpy as np


def crop_wfc3_img(img_file, xmin, xmax, ymin, ymax, output_filename):
    img =  fits.getdata(img_file, 0)
    cropped_img = img[xmin:xmax, ymin:ymax]
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
    crop_wfc3_img('../../images/astrodrizzle_wfc3_ccd_platescale/63.5_rot/final_drz_sci.fits', 814+17, 1858+17, 2376+17, 2444+17, '../../multispec/ccd_multispec/f336w_63.5_cropped.fits')
    crop_wfc3_img('../../images/astrodrizzle_wfc3_ccd_platescale/64.06_rot/final_drz_sci.fits', 814, 1858, 2376, 2444, '../../multispec/ccd_multispec/f336w_64.06_cropped.fits')
    crop_wfc3_img('../../images/astrodrizzle_wfc3_ccd_platescale/63.6_rot/final_drz_sci.fits', 814, 1858, 2376, 2444, '../../multispec/ccd_multispec/f336w_63.6_cropped.fits')



'''
img = fits.getdata('../../images/astrodrizzle_wfc3_ccd_platescale/63.5_rot/final_drz_sci.fits', 0)
cropped_img = img[814+17:1858+17, 2376 + 17:2444 + 17]
pyplot.figure()
pyplot.imshow(cropped_img, interpolation = 'nearest', cmap = 'bone', vmin = 0, vmax = 200)
pyplot.savefig('../../multispec/ccd_multispec/f336w_63.5_cropped.pdf')
hdu = fits.PrimaryHDU(cropped_img)
hdu.header.set('Instrume', 'WFC3')
hdu.header.set('Proposid', 11360)
hdu.header.set('xlim', '(2393, 2461)')
hdu.header.set('ylim', '(831, 1875)')
hdu.header.set('EXPTIME', fits.getval('../../images/astrodrizzle_wfc3_ccd_platescale/63.5_rot/final_drz_sci.fits','exptime', 0))
hdulist = fits.HDUList([hdu])
hdulist.writeto('../../multispec/ccd_multispec/f336w_63.5_cropped.fits', clobber = True)


img = fits.getdata('../../images/astrodrizzle_wfc3_ccd_platescale/64.06_rot/final_drz_sci.fits', 0)
cropped_img = img[814:1858, 2376:2444]
pyplot.figure()
pyplot.imshow(cropped_img, interpolation = 'nearest', cmap = 'bone', vmin = 0, vmax = 200)
pyplot.savefig('../../multispec/ccd_multispec/f336w_64.06_cropped.pdf')
hdu = fits.PrimaryHDU(cropped_img)
hdu.header.set('Instrume', 'WFC3')
hdu.header.set('Proposid', 11360)
hdu.header.set('xlim', '(2376, 2444)')
hdu.header.set('ylim', '(814, 1858)')
hdu.header.set('EXPTIME', fits.getval('../../images/astrodrizzle_wfc3_ccd_platescale/64.06_rot/final_drz_sci.fits','exptime', 0))
hdulist = fits.HDUList([hdu])
hdulist.writeto('../../multispec/ccd_multispec/f336w_64.06_cropped.fits', clobber = True)

img = fits.getdata('../../images/astrodrizzle_wfc3_ccd_platescale/63.6_rot/final_drz_sci.fits', 0)
cropped_img = img[814:1858, 2376:2444]
pyplot.figure()
pyplot.imshow(cropped_img, interpolation = 'nearest', cmap = 'bone', vmin = 0, vmax = 200)
pyplot.savefig('../../multispec/ccd_multispec/f336w_63.6_cropped.pdf')
hdu = fits.PrimaryHDU(cropped_img)
hdu.header.set('Instrume', 'WFC3')
hdu.header.set('Proposid', 11360)
hdu.header.set('xlim', '(2376, 2444)')
hdu.header.set('ylim', '(814, 1858)')
hdu.header.set('EXPTIME', fits.getval('../../images/astrodrizzle_wfc3_ccd_platescale/63.6_rot/final_drz_sci.fits','exptime', 0))
hdulist = fits.HDUList([hdu])
hdulist.writeto('../../multispec/ccd_multispec/f336w_63.6_cropped.fits', clobber = True)
'''