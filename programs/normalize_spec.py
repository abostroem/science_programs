import scipy
from scipy import signal
from scipy.signal import medfilt

from matplotlib import pyplot
import numpy as np
import pdb
from astropy import constants
import numpy as np
import pyfits

def normalize_spectrum(wl, net, kernel_size = 201):
    fig = pyplot.figure(figsize = [22, 7])
    ax = fig.add_subplot(1,1,1)
    wl = wl.flatten()
    net = net.flatten()
    ax.plot(wl, net, 'b')
    background_wl = np.empty((0,))
    background_net = np.empty((0,))
    raw_input('Zoom in on spectrum and press enter to continue ')
    select_another_region = 'y'
    print 'Select background regions '
    while select_another_region != 'n':
        pt1, pt2 = pyplot.ginput(n = 2, timeout = -1)
        x1, y1 = pt1
        x2, y2 = pt2
        background_indx = np.where((wl < max(x1, x2)) & (wl > min(x1, x2)))
        background_wl = np.append(background_wl, wl[background_indx])
        background_net = np.append(background_net, net[background_indx])
        pyplot.plot([wl[background_indx[0][0]], wl[background_indx[0][-1]]], [net[background_indx[0][0]], net[background_indx[0][-1]]], 'm--|', lw = 3)

        select_another_region = raw_input('Would you like to select another background regions? (y), n ')
    background_wl = background_wl.flatten()
    background_net = background_net.flatten()
    sort_indx = np.argsort(background_wl)
    background_wl = background_wl[sort_indx]
    background_net = background_net[sort_indx]
    poly_coeff = np.polyfit(background_wl, background_net, 3)
    continuum = np.polyval(poly_coeff, wl)
    ax.plot(wl, continuum, 'c', lw = 2)
    pyplot.draw()
    return wl, net/continuum, continuum

def apply_redshift(line_dict, v):
    #velocity should be in m/s
    redshift_list = []
    for ilist in line_dict['rest_wavelength']:
        elem_list = []
        for iline in ilist:
            elem_list.append(iline* (1.0 + v/constants.c.value))
        redshift_list.append(elem_list)
    line_dict['redshift_wavelength'] = redshift_list
    return line_dict

def mark_spectrum(line_dict, ax, c, offset = 0.0):
    for label, lines in zip(line_dict['name'], line_dict['redshift_wavelength']):
        for value in lines:
            ax.axvline(value, color = c, linestyle = '--')
        ymin, ymax = ax.get_ylim()
        t = ax.text(lines[0]+offset, ymax - (ymax - ymin)/15.0, label, horizontalalignment = 'right',  rotation = 'vertical',  color = c)

        #t.set_fontweight('bold')
        #t.set_backgroundcolor('w')
        #t.set_rotation('vertical')
    return ax
        
def manually_remove_CR(wl, spec, flt = None):
    fig = pyplot.figure(figsize = [22, 14])
    ax = fig.add_subplot(2,1,1)
    if flt:
        ax2 = fig.add_subplot(2, 1, 2)
        img = pyfits.getdata(flt, 1)
        ax2.imshow(img, interpolation = 'nearest', aspect = 'auto', cmap = 'bone', vmin = 0, vmax = 1000)
    ax.plot(wl, spec)
    pyplot.draw()
    remove_CR_flag = raw_input('Would you like to remove any cosmic-rays? (y), n ')
    while remove_CR_flag != 'n':
        raw_input('Zoom in on cosmic ray and press ENTER when finished ')
        print 'Click on cosmic ray on plot'
        x, y = pyplot.ginput(n = 1, timeout = -1)[0]
        rm_indx = np.where(abs(wl - x) == np.min(abs(wl - x)))[0]
        spec[rm_indx] = (spec[rm_indx - 1] + spec[rm_indx + 1])/2.0
        print 'CR removed at %f' %(wl[rm_indx])
        ax.cla()
        ax.plot(wl, spec)
        pyplot.draw()
        remove_CR_flag = raw_input('Would you like to remove any cosmic-rays? (y), n ')
    pyplot.close()
    return spec





NII_dict = {'name':['N II 3995', 'N II 4041/44'],\
            'rest_wavelength':[[3995],  [4041, 4044]]}

NIII_dict = {'name':['N III 4097','N III 4379' ,'N III 4511/15','N III 4535/414/42','N III 4905'],\
            'rest_wavelength':[[4097],[4379],[4511, 4515],[4634, 4641, 4642],[4905]]}
NIV_dict = {'name':['N IV 4058'], \
            'rest_wavelength': [[4058]]}

NV_dict = {'name':['N V 1238/42','N V 4604/20'],\
            'rest_wavelength':[[1238, 1242],[4604, 4620]]}

HeI_dict = {'name':['He I 4009' ,'He I 4121','He I 4144','He I 4388', 'He I 4471','He I 4713','He I 4921'],\
            'rest_wavelength':[[4009], [4121],[4144],  [4388],[4471],[4713],[4921]]}

HeI_II_dict = {'name':['He I+II 4026'],\
            'rest_wavelength':[[4026]]}

HeII_dict = {'name':['He II 1640', 'He II 4200','He II 4542','He II 4686'],\
            'rest_wavelength':[ [1640], [4200],[4542], [4686]]}

NiIV_dict = {'name':['Ni IV 4058'],\
            'rest_wavelength':[[4058]]} 

CIV_dict = {'name':['C IV 1548/51'], \
            'rest_wavelength':[[1548, 1551]]}

CIII_dict = {'name':['C III 1176','C III 4068/69/70','C III 4187','C III 4647/50/51'],\
            'rest_wavelength':[ [1176],[4068, 4069, 4070],[4187],[4647, 4650, 4651]]}

CII_dict = {'name':['C II 4267'],\
            'rest_wavelength':[[4267]]}

OII_dict = {'name':['O II 4070/72/76','O II 4254','O II 4276/85','O II 4317/20', 'O II 4349','O II 4367','O II 4415/17','O II 4591/4596' ,'O II 4662/76' ,'O II 4699/705'],\
            'rest_wavelength':[[4070, 4072, 4076],[4254],[4276, 4285], [4317, 4320],[4349], [4367], [4415, 4417],[4591, 4596],[4662, 4676],[4699, 4705]]}

OV_dict = {'name':['O V 1371'],\
            'rest_wavelength':[[1371]]}

SiIV_dict = {'name':['Si IV 1393/402','Si IV 4089','Si IV 4116' ,'Si IV 4631'],\
            'rest_wavelength':[[1391, 1402], [4089], [4116], [4631]]}

SiIII_dict = {'name':['Si III 4553/68/75'],\
            'rest_wavelength':[[4553, 4568, 4575]]}

MgII_dict = {'name':['Mg II 4481'],\
            'rest_wavelength':[[4481]]}

SIV_dict = {'name':['S IV 4486','S IV 4504'],\
            'rest_wavelength':[[4486],[4504]]}

balmer_dict = {'name':['H - $epsilon$ 3970', 'H-$\delta$ 4102','H-$\gamma$ 4340','H-$beta$ 4861'],\
            'rest_wavelength':[[3970], [4102], [4340],[4861]]}

dib_dict = {'name':['DIB is 4429','DIB is 4502','DIB is 4727','DIB is 4762','DIB is 4780','DIB is 4881'],\
            'rest_wavelength':[[4429],  [4502], [4727],[4762], [4780],[4881]]}

CaII_dict = {'name': ['Ca II H 3968', 'Ca II K 3933'], \
            'rest_wavelength': [[3968], [3933]]}