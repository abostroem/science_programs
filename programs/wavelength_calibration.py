import pyfits
import numpy as np
import math
import sys
import os
from matplotlib import pyplot
import pdb
import glob
from scipy import constants

oref = '/Users/bostroem/science/oref/'
os.environ['oref'] = '/Users/bostroem/science/oref/'


def get_spec_parameters(spec_file):
    '''
    This function reads in a spectrum and the primary and first extension headers
    Input:
        spec_file: the name of the file
    Output:
        spec: binary fits table
        hdr0: primary header
        hdr1: first extension header
    '''
    spec = pyfits.getdata(spec_file, 1)
    hdr0 = pyfits.getheader(spec_file, 0)
    hdr1 = pyfits.getheader(spec_file, 1)
    return spec,hdr0, hdr1

def get_disp_coeff(hdr, a2center):
    '''
    This function reads in the initial dispersion coefficients and reference aperture from
    the DISPTAB.
    Inputs:
        hdr: Primary header of the spectrum
        a2center: the y location of the extracted spectra
    Outputs:
        disp_coeff: dispersion coefficients
        ref_aper: reference aperture
    '''
    disptab = os.path.join(oref, hdr['disptab'].split('$')[-1])
    tbdata = pyfits.getdata(disptab, 1)
    if a2center == 513.0:
        print 'WARNING: a2center is set to 513.0. It should only have this value if you don\'t know where the spectrum was extracted from'
    row_indx = np.where((tbdata['opt_elem'] == hdr['opt_elem']) & (tbdata['cenwave'] == hdr['cenwave']))
	#a2center should be 1 indexed
    row_indx_a2center = row_indx[0][np.argmin(np.abs(tbdata['a2center'][row_indx] - a2center))]
    disp_coeff = tbdata['coeff'][row_indx_a2center][0:tbdata['ncoeff'][row_indx_a2center]]
    ref_aper = tbdata['ref_aper'][row_indx_a2center]
    return disp_coeff, ref_aper

def get_inangtab_coeff(hdr):
    '''
    Get the correct coefficient for the incidence angle correction. As of the writing of
    this code, only one of the coefficients is non-zero. If this ever changes, this section
    should be rewritten
    
    Input:
        hdr: primary header from spectrum
    Output:
        coeff1: indicent angle correction coefficient
    '''
    inangtab = os.path.join(oref, hdr['inangtab'].split('$')[-1])
    tbdata = pyfits.getdata(inangtab, 1)
    row_indx = np.where(tbdata['opt_elem'] == hdr['opt_elem']) 
    if tbdata['ncoeff2'][row_indx] != 0:
        print 'NCOEFF2 != 0, REWRITE CODE'
        sys.exit()
    coeff1 = tbdata['coeff1'][row_indx][0][0:tbdata['ncoeff1'][row_indx]]
    #pdb.set_trace()
    return coeff1


def get_aperture_offset(hdr, ref_aper):
    '''
    Calculate the correction for using an aperture which is not the reference aperture in
    the DISPTAB.
    Input:
        hdr: primary header from the spectrum file
        ref_aper: the reference aperture read from the DISPTAB
    Output:
        diff_aper_offset1 = aperture offset
    '''
    apdestab = os.path.join(oref, hdr['apdestab'].split('$')[-1])
    tbdata = pyfits.getdata(apdestab, 1)
    row_indx = np.where(tbdata['aperture'] == hdr['aperture'])
    aper_offset1 = tbdata['offset1'][row_indx]
    ref_row_indx = np.where(tbdata['aperture'] == ref_aper)
    ref_aper_offset1 = tbdata['offset1'][ref_row_indx]
    diff_aper_offset1 = aper_offset1 - ref_aper_offset1
    return diff_aper_offset1

    
def modify_disp_coeff(A, c1, s, hdr0, hdr1):  
    '''
    This function corrects the dispersion coefficients as described in Table 3 from 
    STIS ISR 1999-03. Not the MAMA offset is no longer used
    
    Inputs:
        A: dispersion coefficints
        c1: incident angle correction coefficint
        s: aperture offset correction
        hdr0: primary header from the spectrum file
        hdr1: the first extension header from the spectrum file
    output:
        A: corrected dispersion coefficients
    '''

    while len(A) < 7:
        A = np.append(A, 0)
    while len(c1) < 7:
        c1 = np.append(c1, 0)
        
    shifta1 = hdr1['shifta1']

    A = A + c1*s

    A[0] = A[0] + shifta1

    return A

def disp_solution_function(A, pix):
    '''
    Find the wavelength corresponding to each input pixel using the equation 2 from STIS 
    ISR 1999 - 03
    pix = (A[0] + A[3]) + (A[1] + A[5] + A[4])*l  +  (A[2] + A[6])* l**2 
    Normally this is solved using Newton-Raphson. However, since the spectral order is 
    always 1, I can create a quadratic equation and explicitly solve for the roots.
    
    Inputs:
        A: dispersion coefficients
        pix: x pixel numbers
    Outputs:
        wave: wavelengths corresponding to the input pixels
    '''
    c = A[0] + A[3] - pix
    b = A[1] + A[5] + A[4]
    a = A[2] + A[6]
    #pdb.set_trace()
    l_plus = (-b + math.sqrt(b**2 - 4.0*a*c))/(2.0 * a)
    l_minus = (-b - math.sqrt(b**2 - 4.0*a*c))/(2.0 * a)
    if l_plus > 0:
        return l_plus
    l_minus = (-b - math.sqrt(b**2 - 4.0*a*c))/(2.0 * a)
    if l_minus > 0:
        return l_minus
    else:
        print 'All roots less than 0, something is wrong with your code'
        sys.exit()

def helcorr(wave, hdr1):
    '''
    Correct for the motion of Earth around the sun.
    Inputs:
        wave: wavelengths
        hdr1: first extension header from spectrum file
    Outputs:
        wave: corrected wavelengths
    '''
    #Radial velocity is positive when the Earth is receding from the target.
    #--> positive radial velocity is a redshift, negative is a blue shift
    #--> correction must be blueshift for positive radial velocity and red shift for negative radial velocity
    wave = wave* (1.0 - (hdr1['v_helio']*1000.0/constants.c))
    return wave
     
########################################################
#This is the main code of this module which calls to all the other functions
########################################################
                         
def calibrate(composite_file, a2center = 513.0):
    '''
    This code provides a wavelength calibration for manually extracted data from a flt file
    The input file should be a 1D spectrum with a column called net in the first extension as
    a binary fits table. It should have the primary (0th) and first extension header information 
    from the FLT file in the primary and first extension header of the 1D spectrum file.
    
    This calibration is based on STIS ISR 1999-03 and personal communications with Phil Hodge
    The MAMA Offset correction is now performed in the SHIFT1 calculation and not read from
    the table as described in the TIR
    
    This only works for first order spectra.
    
    Input:
        composite_file: spectrum file as detailed above
        a2center: y pixel location from which the data were extracted
    output:
        wave: wavelength array
        spec['net']: net counts from composite_file
    
    Calls to:
        get_spec_parameters
        get_disp_coeff
        get_inangtab_coeff
        get_aperture_offset
        modify_disp_coeff
        helcorr
        
        
    Example:
        wave_star, net_star = calibrate(ifile)
    '''
    spec, hdr0, hdr1 = get_spec_parameters(composite_file)
    disp_coeff, ref_aper = get_disp_coeff(hdr0, a2center = a2center)
    inang_coeff1 = get_inangtab_coeff(hdr0)
    aper_offset1 = get_aperture_offset(hdr0, ref_aper)
    disp_coeff = modify_disp_coeff(disp_coeff, inang_coeff1, aper_offset1, hdr0, hdr1)
    wave = np.empty((0,))
    for pix in np.arange(len(spec['net'].ravel())) + 1.0:
        wave = np.append(wave, disp_solution_function(disp_coeff, pix))
    wave = helcorr(wave, hdr1)
    return wave, spec['net']
 
