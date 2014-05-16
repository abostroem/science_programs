import glob
from matplotlib import pyplot
import pyfits
import numpy as np
import scipy
from scipy import optimize, stats
import pdb

def blink_files():
    for ifile in flist:
        tbdata = pyfits.getdata(ifile, 0)
        tbdata2 = pyfits.getdata('../ORIG/'+ifile[0:9]+'_raw.fits', int(ifile.replace('.', '_').split('_')[-2]))
        pyplot.figure(1, figsize = [20, 20])
        img1 = pyplot.imshow(tbdata, interpolation = 'nearest')
        img2 = pyplot.imshow(tbdata2, interpolation = 'nearest')
        pyplot.title(ifile)
        img1.set_clim(np.median(tbdata) - 10, np.median(tbdata)+10)
        img2.set_clim(np.median(tbdata2) - 10, np.median(tbdata2)+10)
        pyplot.draw()
        pyplot.xlim(431, 630)
        pyplot.ylim(431, 630)
        for i in range(3):
            img2.set_visible(False)
            pyplot.draw()
            time.sleep(1)
            img2.set_visible(True)
            pyplot.draw()
            time.sleep(1)
        time.sleep(3)
        pyplot.close()


#Make and fit histogram
def fit_func(x, a0, a1, a2):
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2)
    return y

def find_mean_stdev(tbdata):
    #pyplot.figure()
    hist_y, hist_x_edge = np.histogram(tbdata.ravel(), bins = np.arange(100)*3+np.mean(tbdata) - 150.0)
    hist_x = (hist_x_edge[1:] - hist_x_edge[:-1])/2.0 + hist_x_edge[:-1]
    p0 = [np.max(hist_y), np.mean(hist_x), 50.0]
    popt, pcov = scipy.optimize.curve_fit(fit_func, hist_x, hist_y, p0 = p0)
    fit = fit_func(hist_x, *popt)
    ##pyplot.plot(hist_x, fit)
    ##pyplot.close()
    return popt[1], popt[2]
    
def create_pdfs(save = True, interactive = False, glob_str = '*rwc*.fits'):
    print glob_str
    flist = glob.glob(glob_str)
    print flist
    mean_before_arr = []
    mean_after_arr = []
    stdev_before_arr = []
    stdev_after_arr = []
    #Make 3 plots
    if interactive is False:
        pyplot.ioff()
    pyplot.gray()
    for ifile in flist:
        tbdata = pyfits.getdata(ifile, 0)
        tbdata2 = pyfits.getdata(ifile.replace('PSUB', 'ORIG').split('_')[0]+'_raw.fits', int(ifile.replace('.', '_').split('_')[-2]))
        mean_before, stdev_before = find_mean_stdev(tbdata2)
        mean_after, stdev_after = find_mean_stdev(tbdata)
        mean_before_arr.append(mean_before)
        mean_after_arr.append(mean_after)
        stdev_before_arr.append(stdev_before)
        stdev_after_arr.append(stdev_after)
        if mean_before < 3000.0:
            fig1 = pyplot.figure(1, figsize = [20, 10])
            ax1 = fig1.add_subplot(1, 3, 1)
            ax2 = fig1.add_subplot(1, 3, 2)
            ax3 = fig1.add_subplot(1, 3, 3)
            img1 = ax2.imshow(tbdata, interpolation = 'nearest')
            img2 = ax1.imshow(tbdata2, interpolation = 'nearest')
            img3 = ax3.imshow(tbdata2 - tbdata, interpolation = 'nearest')
    
            ax1.set_xlabel('Before  %5.2f $\pm$ %4.2f' %(mean_before, stdev_before))
            ax2.set_xlabel('After %5.2f $\pm$ %4.2f' %(mean_after, stdev_after))
            ax3.set_xlabel('Diff %5.2f $\pm$ %4.2f' %(np.mean(tbdata2-tbdata), scipy.stats.tstd(tbdata2-tbdata)))
            ax2.set_title(ifile)
            img1.set_clim(np.median(tbdata) - 0.02*np.median(tbdata), np.median(tbdata)+0.02*np.median(tbdata))
            img2.set_clim(np.median(tbdata2) - 0.02*np.median(tbdata), np.median(tbdata2)+0.02*np.median(tbdata))
            pyplot.draw()
            ax1.set_xlim(431, 630)
            ax1.set_ylim(431, 630)
            ax2.set_xlim(431, 630)
            ax2.set_ylim(431, 630)
            ax3.set_xlim(431, 630)
            ax3.set_ylim(431, 630)
            if interactive is True:
                pyplot.draw()
                raw_input('Examine plots and press enter to continue')
            not_reviewed = False
            if save is True:
                pyplot.savefig('../PLOTS/%s' %(ifile.replace('.fits', '.pdf')))
            pyplot.close()
        else:
            print 'mean > 3000 and removing the herringbone pattern will not affect the overall SNR for %s' %(ifile)
            not_reviewed = True
    if interactive is False:
        pyplot.ion()
    return not_reviewed
    
    


if __name__ == "__main__":
    create_pdfs()
    blink_files()