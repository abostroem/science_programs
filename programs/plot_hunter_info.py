from astropy.io import ascii
import pyfits
from matplotlib import pyplot
import numpy as np
import pdb

class MyFigure:
    def __init__(self):
        self.fig = pyplot.figure(1, figsize = [25, 20])
        '''
        self.ax1 = self.fig.add_subplot(1, 4, 1)
        self.ax2 = self.fig.add_subplot(1, 4, 2)
        self.ax3 = self.fig.add_subplot(1, 4, 3)
        self.ax4 = self.fig.add_subplot(1, 4, 4)  
        '''
        self.ax1 = self.fig.add_subplot(1, 3, 1)
        self.ax2 = self.fig.add_subplot(1, 3, 2)
        self.ax3 = self.fig.add_subplot(1, 3, 3)
        #'''
def display_img(img):

    fig_obj = MyFigure()

    img = pyfits.getdata(img, 1)
    im1 = fig_obj.ax1.imshow(img, cmap = 'bone', vmin = 0, vmax = 1000)
    im2 = fig_obj.ax2.imshow(img, cmap = 'bone', vmin = 0, vmax = 1000)
    im3 = fig_obj.ax3.imshow(img, cmap = 'bone', vmin = 0, vmax = 4000)
    im4 = fig_obj.ax4.imshow(img, cmap = 'bone', vmin = 0, vmax = 1000)

    fig_obj.ax1.set_xlim(340, 420)
    fig_obj.ax2.set_xlim(340, 420)
    fig_obj.ax3.set_xlim(340, 420)
    fig_obj.ax4.set_xlim(340, 420)
    
    fig_obj.ax1.set_ylim(575,700)
    fig_obj.ax2.set_ylim(450,575)
    fig_obj.ax3.set_ylim(325,450)
    fig_obj.ax4.set_ylim(200, 325)
    

    return fig_obj

def display_img2(img):

    fig_obj = MyFigure()

    img = pyfits.getdata(img, 1)
    im1 = fig_obj.ax1.imshow(img, cmap = 'Blues_r', vmin = 0, vmax = 1000)
    im2 = fig_obj.ax2.imshow(img, cmap = 'Blues_r', vmin = 0, vmax =4000)
    im3 = fig_obj.ax3.imshow(img, cmap = 'Blues_r', vmin = 0, vmax = 1000)

    fig_obj.ax1.set_xlim(340, 420)
    fig_obj.ax2.set_xlim(340, 420)
    fig_obj.ax3.set_xlim(340, 420)

    fig_obj.ax1.set_ylim(533,700)
    fig_obj.ax2.set_ylim(367,533)
    fig_obj.ax3.set_ylim(200,367)
    

    return fig_obj

def label_stars(filename, column, fig_obj):
    data_table = ascii.read(filename, names = ['number', 'x', 'y', 'f336w', 'f336w_err','f555w', 'f555w_err', 'f814', 'f814_err', 'color1', 'color1_err', 'color2', 'color2_err' ], data_start = 0, data_end = 1000)
    
    #ax_list = [fig_obj.ax1, fig_obj.ax2, fig_obj.ax3, fig_obj.ax4]
    #ax_list = [fig_obj.ax1, fig_obj.ax2]
    ax_list = [fig_obj.ax1, fig_obj.ax2, fig_obj.ax3]
    for ax in ax_list:
        ax_xlims = ax.get_xlim()
        ax_ylims =ax.get_ylim()
        indx = np.where((data_table['x'] < ax_xlims[1]) & (data_table['x'] > ax_xlims[0]) & (data_table['y'] > ax_ylims[0]) & (data_table['y'] < ax_ylims[1]) & (data_table['f336w'] < 100.0))
        ax.plot(data_table['x'][indx]-1, data_table['y'][indx]-1, '.', color = 'r')
        for x, y, s in zip(data_table['x'][indx]-.9, data_table['y'][indx]-.9, data_table[column][indx]):
            ax.text(x, y, s, color = 'r')
    fig_obj.fig.suptitle('WFPC2 image with Hunter stars (1st 1000). Labels = %s' %(column))
    

def label_50_brighest_stars(filename, column, fig_obj):
    data_table_all = ascii.read(filename, names = ['number', 'x', 'y', 'f336w', 'f336w_err','f555w', 'f555w_err', 'f814', 'f814_err', 'color1', 'color1_err', 'color2', 'color2_err' ], data_start = 0)

    data_table_sort = data_table_all[np.argsort(data_table_all['f555w'])]
    indx_all = np.where((data_table_sort['x'] > 340) & (data_table_sort['x'] <= 420) & (data_table_sort['y'] > 200) & (data_table_sort['y'] < 700) & (data_table_sort['f555w'] < 100))[0]
    data_table = data_table_sort[indx_all][0:51]
    ax_list = [fig_obj.ax1, fig_obj.ax2, fig_obj.ax3]
    rhodes_list = [9, 36, 19, 62, 50, 48, 55, 47, 40, 35, 24, 58]
    aas_list = [5, 66, 32, 46, 45, 52, 92, 73, 68, 49, 86, 39]
    counter = 0
    for ax in ax_list:
        ax.plot(data_table_sort['x'][indx_all][0:1000]-1, data_table_sort['y'][indx_all][0:1000]-1, '.', color = 'y')
        print  data_table_sort['f555w'][indx_all][0:1000]
        ax_xlims = ax.get_xlim()
        ax_ylims =ax.get_ylim()
        indx = np.where((data_table['x'] < ax_xlims[1]) & (data_table['x'] > ax_xlims[0]) & (data_table['y'] > ax_ylims[0]) & (data_table['y'] < ax_ylims[1]) & (data_table['f336w'] < 100.0))
        ax.plot(data_table['x'][indx]-1, data_table['y'][indx]-1, '.', color = 'r')
        for x, y, s1, s2 in zip(data_table['x'][indx]-.9, data_table['y'][indx]-.9, data_table['f555w'][indx], data_table['number']):
                ax.text(x, y, s1, color = 'r')
                ax.text(x-3, y-2, s2, color = 'DeepPink')
                if s2 in rhodes_list:
                    ax.plot(x, y, 's', color = 'r', markersize =7)
                if s2 in aas_list:
                    ax.plot(x, y, '^', color = 'g', markersize = 10, markeredgecolor = 'g')

        
    pdb.set_trace()
    data_table = data_table_sort[indx_all][51:101]
    for ax in ax_list:
        ax_xlims = ax.get_xlim()
        ax_ylims =ax.get_ylim()
        indx = np.where((data_table['x'] < ax_xlims[1]) & (data_table['x'] > ax_xlims[0]) & (data_table['y'] > ax_ylims[0]) & (data_table['y'] < ax_ylims[1]) & (data_table['f336w'] < 100.0))
        ax.plot(data_table['x'][indx]-1, data_table['y'][indx]-1, '.', color = 'y')
        for x, y, s1, s2 in zip(data_table['x'][indx]-.9, data_table['y'][indx]-.9, data_table['f555w'][indx], data_table['number']):
                ax.text(x, y, s1, color = 'GoldenRod')
                ax.text(x-3, y-2, s2, color = 'DarkOrange')
                if s2 in rhodes_list:
                    ax.plot(x, y, 's', color = 'Olive', markersize = 5)
            
            
    #fig_obj.fig.suptitle('WFPC2 image with Hunter stars (1st 1000). Labels = %s' %(column))
    pdb.set_trace()



if __name__ == "__main__":
    fig_obj = display_img2('..//multispec/multi_spec_files/r136_f555w_wfpc2_images/align_w_wfc3_wcs/u25y0107t_c0m.fits')
    label_50_brighest_stars('../multispec/multi_spec_files/hunter_ApJ1995_488_179.dat', 'f555w', fig_obj)
    pdb.set_trace()
    label_stars('../multispec/multi_spec_files/hunter_ApJ1995_488_179.dat', 'number', fig_obj)
    pdb.set_trace()
    pyplot.close()
    fig_obj = display_img2('..//multispec/multi_spec_files/r136_f555w_wfpc2_images/align_w_wfc3_wcs/u25y0107t_c0m.fits')
    label_stars('../multispec/multi_spec_files/hunter_ApJ1995_488_179.dat', 'f555w', fig_obj)
    pdb.set_trace()
