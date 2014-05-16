import pyfits
import numpy as np
from matplotlib import pyplot
import os

def make_plot(filename):
    img = pyfits.getdata(filename, 1)
    spec = np.sum(img[:, 400:600], axis = 1)
    pyplot.plot(spec, np.arange(np.shape(spec)[0]))
    pyplot.ylim(450, 580)
    pyplot.savefig('/user/bostroem/science/12465_otfr20121109/ccd/cross_disp_profiles/%s' %(filename.replace('fits', 'pdf')))
    pyplot.close()



if __name__ == "__main__":
    os.chdir('/user/bostroem/science/12465_otfr20121109/ccd/')
    make_plot('NW1_3936_combined_img.fits')
    make_plot('NW2_3936_combined_img.fits')
    make_plot('NW3_3936_combined_img.fits')
    make_plot('NW4_3936_combined_img.fits')
    make_plot('NW5_3936_combined_img.fits')
    make_plot('NW6_3936_combined_img.fits')
    make_plot('NW7_3936_combined_img.fits')
    make_plot('NW8_3936_combined_img.fits')
    make_plot('SE1_3936_combined_img.fits')
    make_plot('SE2_3936_combined_img.fits')
    make_plot('SE3_3936_combined_img.fits')
    make_plot('SE4_3936_combined_img.fits')
    make_plot('SE5_3936_combined_img.fits')
    make_plot('SE6_3936_combined_img.fits')
    make_plot('SE7_3936_combined_img.fits')
    make_plot('SE8_3936_combined_img.fits')
    make_plot('SE9_3936_combined_img.fits')
