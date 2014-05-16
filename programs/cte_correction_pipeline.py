import os
import shutil
import glob
import sys
try:
    import pyfits
except:
    from astropy.io import fits as pyfits
import numpy as np
import pdb
import time
import datetime
from astropy.time import Time
import stistools

#Set paths
#Eventually this can be done with argparse
#set cur_dir, refstis_dir, 
#What do I use programs for?

cur_dir = '/Users/bostroem/science/cte/2012_04'
refstis_dir = '/Users/bostroem/science/programs/refstis/'
ctecorr_dir = '/Users/bostroem/science/programs/ctecorr/stis_mac/lib/stistools'
os.environ['myref'] = os.path.join(cur_dir, 'reffiles')+'/'

os.environ['oref'] = '/Users/bostroem/science/oref/'

herringbone_dir = os.path.join(cur_dir, 'herringbone_correction')
sys.path.append('/Users/bostroem/science/programs')
sys.path.append(herringbone_dir)
sys.path.append(refstis_dir)
sys.path.append(ctecorr_dir)
pctetab = 'stis_apr_2012_pcte.fits'
import REFSTIS_wrapper
import batch_run_autofilet
from REFSTIS_functions import divide_anneal_month
import StisPixCteCorr





def apply_herringbone_correction(input_dir, input_flist, output_dir = None):
    '''
    This function will copy files from an input folder to the ORIG folder in the herringbone
    correction folder. It will then run the herringbone correction script in batch mode on those files
    and copy the final result from the CLEAN folder back to the original folder
    '''
    if not output_dir:
        output_dir = input_dir
    #Clean out herringbone directory
    print 'Removing old files from herringbone correction folder'
    clean_herringbone_folders()
    #Copy files to the ORIG directory in the herringbone_correction folder
    for ifile in input_flist:
        shutil.copy(os.path.join(input_dir, ifile), 'herringbone_correction/ORIG/%s' %(ifile))
    #Run autofilet
    os.chdir(herringbone_dir)
    print 'Starting Herringbone Correction'
    batch_run_autofilet.autofilet()
    batch_run_autofilet.rebuild_files(interactive = False)
    #Move original raw files to the orig folder and corrected files to input_dir
    os.chdir(os.path.join(cur_dir, input_dir))
    if not os.path.exists('orig'):
        os.mkdir('orig')
    for ifile in input_flist:
        shutil.copy(ifile, os.path.join(cur_dir, input_dir, 'orig/%s' %(ifile))) #copy uncorrected to orig
        shutil.copy(os.path.join(herringbone_dir, 'CLEAN/%s' %(ifile)), os.path.join(cur_dir, output_dir, ifile)) #copy corrected to input_dir
    os.chdir(cur_dir)
    print 'HERRINGBONE CORRECTION APPLIED '
    print '\tinput directory = %s' %(input_dir)
    print '\toutput directory = %s' %(output_dir)

def clean_herringbone_folders():
    cur_dir = os.getcwd()
    os.chdir(herringbone_dir)
    flist = glob.glob('PSUB/*')+glob.glob('CLEAN/*') + glob.glob('TABLES/*') + glob.glob('ORIG/*') + glob.glob('FAILED/*')
    for ifile in flist:
        os.remove(ifile)    
    os.chdir(cur_dir)

def create_reference_file(input_dir, input_flist):
    '''
    This function will point Justin's refstis scripts to a folder containing darks or biases
    and have it create a superdark or super bias. The script will automatically try to create darks
    and biases for with all files in a directory. Biases and darks should be kept in separate 
    folders. This script automatically divided observations into individual weeks
    '''
    date = []
    for ifile in input_flist:
        date.append(pyfits.getval(os.path.join(input_dir, ifile), 'texpstrt', 0))
        date.append(pyfits.getval(os.path.join(input_dir, ifile),'texpend', 0))
    filetype = pyfits.getval(os.path.join(input_dir, ifile), 'targname', 0)
    date = np.array(date)
    month_begin = np.min(date)
    month_end = np.max(date)
    if filetype == 'DARK':
        REFSTIS_wrapper.separate_obs(input_dir, month_begin, month_end, input_filetype = 'flc', move_files = False)
    else:
        REFSTIS_wrapper.separate_obs(input_dir, month_begin, month_end, move_files = False)
    make_reffiles(input_dir, filetype)

    #Copy reference file from each week folder to the reffiles folder
    os.path.walk(input_dir, copy_reference_files_to_reffiles, '')   

def make_reffiles(root_folder, filetype):
    """ Make all refrence files for a given folder

    This functions is very specific to the REFSTIS pipeline, and requires files
    and folders to have certain naming conventions.  

    """
    print '#-----------------------------#'
    print '#  Making all ref files for   #'
    print  root_folder
    print '#-----------------------------#'

    if not os.path.exists( root_folder ): raise IOError( 'Root folder: %s does not exist' %(root_folder))
    gain_folders, week_folders = REFSTIS_wrapper.pull_out_subfolders( root_folder )
    if filetype.upper() == 'BIAS':
        REFSTIS_wrapper.refbias(week_folders)
        #Basebias is not currently supported   
        #REFSTIS_wrapper.basebias(root_folder, gain_folders)
    elif filetype.upper() == 'DARK':
        REFSTIS_wrapper.basedark(root_folder, input_filetype = 'flc', set_biasfile = False)
        REFSTIS_wrapper.weekdark(week_folders, input_filetype = 'flc', set_biasfile = False)
    else:
        sys.exit('ERROR: File type not recognized')

def copy_reference_files_to_reffiles(junk, dirname, junk2):
    '''
    This function copies super bias and super dark reference files to the reffiles folder
    '''
    flist = glob.glob(os.path.join(dirname, '*bia.fits')) + glob.glob(os.path.join(dirname, '*drk.fits'))
    if len(flist) > 0:
        for ifile in flist:
            shutil.copy(ifile, os.path.join(cur_dir, 'reffiles', ifile.split('/')[-1]))
            print 'Copying %s to reffiles folder' %(ifile)

def calibrate_data_with_cte_correction(input_dir, input_flist):
    '''
    This function will:
    1. copy relevant reference files
    2. update relevant header keywords
    3. Start running calstis with BIASCORR, BLEVCORR, and DQICORR
    4. Run the CTE correction code
    5. Finish running calstis on the CTE corrected files with CTECORR turned off
    '''
    #copy_relevant_reference_files(reffile_list, output_dir)
    update_reffile_keywords(input_dir, input_flist, 'bia')
    try:
        update_reffile_keyword(input_dir, input_flist, 'drk')
    except:
        print 'WARNING: Not updating DARKFILE for files in %s' %(input_dir)
    add_pctetab_to_headers(input_dir, input_flist)
    run_calstis_part1(input_dir, input_flist)
    run_cte_correction_code(input_dir, input_flist)
    run_calstis_part2(input_dir, input_flist)


def update_reffile_keywords(input_dir, input_list, filetype):
    '''
    This function takes the files in input_dir/input_flist and updates
    the reference files associated with filetype 
    Inputs:
        input_dir: directory of files whose headers you want to update
        input_list: names of files whose headers you want to update
        filetype: either drk or bia - the reference file whose name you want to update
    Outputs:
        modifies files in input_flist
    '''
    #Get dates for reference file application as a function of observation mode
    mode_dict = determine_correct_reference_files(input_dir, input_list, filetype)
    #Determine which keyword to update
    assert (filetype is 'bia') or (filetype is 'drk'), 'ERROR: Unknown filetype %s, please use bia or drk' %(filetype)
    filetype_dict = {'bia':'BIASFILE', 'drk':'DARKFILE'}
    keyword = filetype_dict[filetype]
    #Update headers
    for ifile in input_list:
        hdr0 = pyfits.getheader(os.path.join(input_dir, ifile), 0)
        gain = hdr0['ccdgain']
        binaxis1 = hdr0['binaxis1']
        binaxis2 = hdr0['binaxis2']
        expstart = hdr0['texpstrt']
        useafter_dates = mode_dict[gain][binaxis1][binaxis2].keys()
        useafter_dates.sort()
        useafter_dates = np.array(useafter_dates)
        date_diff = expstart - useafter_dates 
        date_indx = np.where((date_diff > 0))
        reffile = mode_dict[gain][binaxis1][binaxis2][useafter_dates[date_indx[0][-1]]]
        pyfits.setval(os.path.join(input_dir, ifile), keyword, value = reffile, ext = 0)

def determine_correct_reference_files(input_dir, input_list, filetype):
    '''
    This function creates a nested dictionary which assigns a reference file name (from 
    the reffiles folder) to each gain, binaxis1, binaxis2, useafter combination
    Inputs:
        input_dir: directory of files whose headers you want to update
        input_list: names of files whose headers you want to update
        filetype: either drk or bia - the reference file whose name you want to update
    Output:
        mode_dict: a nested dictionary which hich assigns a reference file name (from 
    the reffiles folder) to each gain, binaxis1, binaxis2, useafter combination

    This code assumes that you have set the environment variable myref to point to the 
    reffiles folder
    '''
    dates = []
    for ifile in input_list:
        filename = os.path.join(input_dir, ifile)
        dates.append(pyfits.getval(filename, 'texpstrt', 0))
        dates.append(pyfits.getval(filename, 'texpend', 0))
    data_start = min(dates)
    data_end = max(dates)
    anneal_weeks_4 = divide_anneal_month(data_start, data_end, '/Users/bostroem/science/cte/', 4)
    anneal_weeks_2 = divide_anneal_month(data_start, data_end, '/Users/bostroem/science/cte/', 2)
    #nested dictionary: gain, binaxis1, binaxis2, week
                #gain
    mode_dict = {1:
                    #binaxis1
                    {1:
                        #binaxis2
                        #1x1 gain = 1
                        {1:{round(anneal_weeks_4[0][0], 4):'', round(anneal_weeks_4[1][0], 4):'', round(anneal_weeks_4[2][0], 4):'', round(anneal_weeks_4[3][0], 4):''},   
                        # 1x2 gain = 1
                        2:{round(anneal_weeks_2[0][0], 4):'', round(anneal_weeks_2[1][0], 4):''}},    
                    2: 
                        #2x1 gain = 1
                        {1:{round(anneal_weeks_2[0][0], 4):'', round(anneal_weeks_2[1][0], 4):''},  
                        #2x2 gain = 1
                        2:{round(anneal_weeks_2[0][0], 4):'', round(anneal_weeks_2[1][0], 4):''}}}, 
                4:  
                    #binaxis1
                    {1:
                        #binaxis2
                        #1x1 gain = 4
                        {1:{round(anneal_weeks_2[0][0], 4):'', round(anneal_weeks_2[1][0], 4):''}}}}   
    reffile_list = glob.glob(os.path.join(os.environ['myref'], '*_%s*.fits' %(filetype)))
    for reffile in reffile_list:
        hdr0 = pyfits.getheader(reffile, 0)
        gain = hdr0['ccdgain']
        binaxis1 = hdr0['binaxis1']
        binaxis2 = hdr0['binaxis2']
        useafter = convert_useafter_to_mjd(hdr0['useafter'])
        mode_dict[gain][binaxis1][binaxis2][round(useafter, 4)] = 'myref$%s' %(reffile.split('/')[-1])
                       
    return mode_dict

def convert_useafter_to_mjd(useafter_date):
    useafter_reformat = datetime.datetime.strptime(useafter_date, '%b %d %Y %X').strftime('%Y-%m-%d %X')
    useafter_mjd = Time(useafter_reformat, format = 'iso', scale = 'utc').mjd
    return useafter_mjd
        

def remove_calibration_products(filename):
    try:
        os.remove(filename.replace('raw', 'flt'))
    except:
        pass
    try:
        os.remove(filename.replace('raw', 'crj'))
    except:
        pass
    try:
        os.remove(filename.replace('raw', 'flc'))
    except:
        pass
def run_calstis_part1(input_dir, input_flist):
    '''
    This function runs BIASCORR, BLEVCORR and DQICORR on the files in input_flist
    '''
    log = open('calstis_log.txt', 'a')
    log.write('calstis_part1 - basic2d')
    log.close()
    for ifile in input_flist:
        filename = os.path.join(input_dir, ifile)
        remove_calibration_products(filename)

        stistools.basic2d.basic2d(filename, dqicorr = 'perform', blevcorr = 'perform', biascorr = 'perform',
                atodcorr = 'omit', doppcorr = 'omit', lorscorr = 'omit', glincorr = 'omit', lflgcorr = 'omit', 
                darkcorr = 'omit', flatcorr = 'omit', shadcorr = 'omit', photcorr = 'omit', trailer = 'calstis_log.txt')

def run_cte_correction_code(input_dir, input_flist):
    '''
    This function will run the CTE correction code StisPixCteCorr.py on each file in 
    input_flist
    '''
    
    for ifile in input_flist:
        filename = os.path.join(input_dir, ifile)
        StisPixCteCorr.CteCorr(filename.replace('raw', 'flt'),outFits = filename.replace('raw', 'flc') )
        pyfits.setval(filename.replace('raw', 'flc'), 'CTECORR', ext = 0, value = 'COMPLETE')

def run_calstis_part2(input_dir, input_flist):
    '''
    This function finishes running calstis on the CTE corrected data. Note: CTECORR should be
    set to OMIT or COMPLETE
    '''
    #biacorr, blevcorr, and dqicorr are set to complete after part1 is called. I can call calstis
    log = open('calstis_log.txt', 'a')
    log.write('calstis_part2 - calstis')
    log.close()
    for ifile in input_flist:
        filename = os.path.join(input_dir, ifile)
        assert pyfits.getval(filename.replace('raw', 'flc'), 'ctecorr', 0).upper() != 'PERFORM', 'CTECORR = PERFORM in %s' %(filename)
        wavecal = os.path.join(input_dir, pyfits.getval(filename, 'wavecal', 0))
        stistools.calstis.calstis(filename.replace('raw', 'flc'), trailer = 'calstis_log.txt', wavecal = wavecal)

def add_pctetab_to_headers(input_dir, input_flist):
    for ifile in input_flist:
        pyfits.setval(os.path.join(input_dir, ifile), 'PCTETAB', ext = 0, value = 'myref$%s' %(pctetab))


def make_biasfile():
    start = time.time()
    #Herringbone correct Bias files
    os.chdir('bias')
    flist_bias = glob.glob('*raw.fits')
    os.chdir(cur_dir)
    apply_herringbone_correction('bias', flist_bias)

    
    #Create Superbias
    create_reference_file('bias', flist_bias)
    
    end = time.time()
    print 'RUNTIME BIAS = %f minutes' %((end - start)/60.0) 

def make_darkfile():
    start = time.time()
    #Herringbone correct Dark files
    os.chdir('dark')
    flist_dark = glob.glob('*raw.fits')
    os.chdir(cur_dir)
    apply_herringbone_correction('dark', flist_dark)
    
    #CTE correct and calibrate Dark files with new superbias
    calibrate_data_with_cte_correction(os.path.join(cur_dir, 'dark'), flist_dark)
    
    
    #Create Superdark
    create_reference_file('dark', flist_dark)
    end = time.time()
    print 'RUNTIME DARK = %f minutes' %((end - start)/60.0)   

def correct_science_data():
    start = time.time()
    #Herringbone correct science data
    os.chdir('sci_data')
    flist_sci = glob.glob('*raw.fits')
    os.chdir(cur_dir)
    apply_herringbone_correction('sci_data', flist_sci)
    #CTE correcte and calibrate science data with new superbias and superdark
    calibrate_data_with_cte_correction(os.path.join(cur_dir, 'sci_data'), flist_sci)
    end = time.time()
    print 'RUNTIME Science= %f minutes' %((end - start)/60.0) 




if __name__ == "__main__":
    os.chdir(cur_dir)
    #make_biasfile()
    #make_darkfile()
    correct_science_data()



  


    
