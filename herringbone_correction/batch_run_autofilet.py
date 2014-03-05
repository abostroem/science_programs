try:
    import pyfits
except:
    from astropy.io import fits as pyfits
import glob
import os
import subprocess
import time
import numpy as np
import pdb
import shutil
import glob
import sys
from copy import deepcopy
import argparse 

import herringbone_correction_confirmation_plots
from herringbone_correction_confirmation_plots import create_pdfs




def get_extn_nums(ifile):
    '''
    Get the numbers of the science extensions in a given file
    Inputs:
        ifile: name of input file to query for the number of extensions
    '''
    nextend = pyfits.getval(ifile, 'nextend', 0)
    all_ext = np.arange(nextend)+1
    sci_ext = all_ext[::3]
    return sci_ext

def write_idl_file(ifile, output_file, table_file, extn, xord, yord):
    '''
    This code writes the IDL file which is calls autofilet.pro with the appropriate input parameters
    Inputs:
        ifile: input raw files name
        output_file: output file name, this file will be put in the PSUB directory
        table_file: table output file, this file will be put in the TABLES directory
        extn: file extension to be corrected
        xord: xord input to autofilet
        yord: yord input to autofilet
    Output:
        file tmp_idl.pro is writen to disk
    '''

    ofile = open('tmp_idl.pro', 'w')
    #ofile.write('.compile autofilet\n')
    ofile.write('.run autofilet\n')
    ofile.write('autofilet, "%s", "%s", "%s", exten = %i, xorder = %i, yorder = %i, nrej = 5 \n' %(ifile, output_file, table_file, extn, xord, yord))
    ofile.write('exit\n')
    ofile.close()

def write_inventory(ifile, inventory, xord, yord, ext_list):
    '''
    Writes the file inventory.txt. This is equivalent to Jansen's _INVENTORY file
    Inputs:
        ifile: name of input fits file
        inventory: open file handler
        xord: xord output from autofilet
        yord: yord output from autofilet
        ext_list: list of extensions modified
    Output:
        writes inventory line to inventory.txt

    '''
    hdr0 = pyfits.getheader(ifile, 0)
    hdr1 = pyfits.getheader(ifile, 1)
    for ext in ext_list:
        inventory.write('{:<12s} {:<22s} {:<7.2f} {:<2d} {:<5s} {:<4s} {:<8s} {:<9s} {:<7d} {:<7d} {:<14s} {:<15s} {:^+{width}f} {:^+{width}f} {:^+{width}f} {:^{width}f} {:^+{width}f} {:^+{width}f} {:^+{width}f} {:^+{width}f} {:^+{width}f} {:^+{width}f} \n'.format( \
            ifile.split('/')[-1]+'[%d]' %(ext), \
            hdr0['tdateobs']+'T'+hdr0['ttimeobs'], \
            hdr0['texptime'], \
            int(hdr0['ccdgain']), \
            str(hdr0['sizaxis1'])+'x'+str(hdr0['sizaxis2']), \
            str(hdr0['binaxis1'])+'x'+str(hdr0['binaxis2']), \
            hdr0['opt_elem'], \
            hdr0['aperture'], \
            int(hdr0['minwave']), \
            int(hdr0['maxwave']), \
            hdr0['obstype'], \
            hdr0['targname'], \
            np.round(hdr0['ra_targ'], decimals = 4), \
            np.round(hdr0['dec_targ'], decimals = 4), \
            np.round(hdr1['ra_aper'], decimals = 4), \
            np.round(hdr1['dec_aper'], decimals = 4), \
            np.round(hdr1['pa_aper'], decimals = 4), \
            np.round(hdr1['orientat'], decimals = 4), \
            np.round(hdr1['sunangle'], decimals = 4), \
            np.round(hdr1['moonangl'], decimals = 4), \
            np.round(hdr1['sun_alt'], decimals = 4), \
            np.round(hdr1['occdhtav'], decimals = 4), width = 12))

    return inventory



def autofilet():
    '''
    This code replaces the first part of Jansen's readme file. It identifies the correct file type,
    sets parameters, and runs the IDL routine, stopping it if it takes more than 3 minutes to process.
    A warning is printed to the log file in this case. The log (autofilet.log) is appended to if it 
    exists. A note is inserted indicating the beginning of a new run. The code will not process a 
    file/ext combination which already exists in the PSUB folder. This assumse that you have directories 
    called ORIG, TABLES, and PSUB in your current directory.

    the file inventory.txt is written or appended in the same way as autofilet.log. This provides 
    information about each exposure in the ORIG folder

    Calling autofilet will run all of the code

    Calls to:
        write_inventory
        write_idl_file
    '''
    #List of original files
    orig_list = glob.glob('ORIG/*raw.fits') + glob.glob('ORIG/*wav.fits')
    #Open Log File, note if this is a new run (log is appended, not rewritten)
    if os.path.exists('autofilet.log'):
        autofilet_log = open('autofilet.log', 'a', 0)
        autofilet_log.write('NEW RUN STARTS HERE\n')
        autofilet_log.flush()
    else:
        autofilet_log = open('autofilet.log', 'w', 0)
    autofilet_log.close
    #Open inventory file
    if os.path.exists('inventory.txt'):
        inventory = open('inventory.txt', 'a')
    else:
        inventory = open('inventory.txt', 'w')
    inventory.write('# image dateobs texp gain size bin optelem aperture lamb1 lamb2 obstype targname ra_targ dec_targ ra_aper dec_aper pa_aper orientat sunangle moonangle sun_alt occdhtav\n')
    inventory.write('#\n')
    inventory.write("#       image                dateobs        t_exp gain    size   bin optelem aperture  lamb1 lamb2 ___obstype___ ________targname________    ra_targ   dec_targ   ra_aper   dec_aper    pa_aper   orientat   sunangle    moonangle    sun_alt   occdhtav\n") 
    inventory.write("# root_typ.fits[extn] [YYYY-MM-DDThh:mm:ss] [sec][e/adu][pixels]                       [\AA] [\AA]                                            [deg]      [deg]     [deg]      [deg]      [deg]      [deg]      [deg]       [deg]       [deg]     [deg C]\n") 
    inventory.write("#========================================================================================================================================================================================================================================================\n")
    #Loop over all raw files
    ()
    for ifile in orig_list[::-1]:
        ext_array = get_extn_nums(ifile) #get the ext numbers of the science extensions; it is assumed that this is every third one starting with 1
        #Set different inputs for different types of files
        if ifile.split('_')[1] == 'wav.fits':
            xord = 6
            yord = 15
        elif pyfits.getval(ifile, 'targname', 0) in ['BIAS', 'DARK']:
            xord = 15
            yord = 4
        elif ('IMAG' in pyfits.getval(ifile, 'obstype', 0)) & ('ACQ' not in pyfits.getval(ifile, 'obstype', 0)):
            xord = 6
            yord = 6
        elif ('SPEC' in pyfits.getval(ifile, 'obstype', 0)) & ('ACQ' not in pyfits.getval(ifile, 'obstype', 0)):
            xord = 15
            yord = 29
        else:
            print 'could not identify appropriate file for %s' %(ifile)
        #write information to inventory file
        inventory = write_inventory(ifile, inventory, xord, yord, ext_array)
        #Loop over science extensions in each file
        orig_xord = xord
        orig_yord = yord
        for extn in ext_array:
            xord = orig_xord
            yord = orig_yord
            success = False       
            nloop = 0
            output_file = ifile.replace('ORIG', 'PSUB').replace('raw', 'rwc_%i' %(extn)).replace('wav', 'wav_rwc_%i'%(extn))
            if not os.path.exists(output_file):
                
                table_file = ifile.replace('ORIG', 'TABLES').replace('raw.fits', 'fps_%i' %(extn)).replace('wav.fits', 'wav_fps_%i'%(extn)) 
                #Call IDL routine autofilet with the appropriate inputs
                while not success:
                    nloop += 1
                    #sys.stdout = autofilet_log
                    write_idl_file(ifile, output_file, table_file, extn, xord, yord)
                    autofilet_call = subprocess.Popen('/grp/software/Linux/itt/idl/idl81/bin/idl -quiet tmp_idl.pro', shell = True) 
                    #sys.stdout = sys.__stdout__                
                    nwait = 0
                    while (autofilet_call.poll() is None) and (nwait < 7):  #if the code hasn't finished in 3 minutes, terminate it 
                        time.sleep(30)
                        nwait += 1
                        if nwait == 7:
                            autofilet_call.terminate()  #check to see if terminate or kill is better here; this ends the while auto... loop
                    
                            #15  4           <-- if a hang occurs, retry with   17  4 ,  17  5 ,  18  5
                            # 6  6           <-- if a hang occurs, retry with    7  7 ,   6  8 ,   7  8
                            #15 29           <-- if a hang occurs, retry with   15 31 ,  17 29 ,  16 30
                            # 6 15           <-- if a hang occurs, retry with    6 17 ,   7 15 ,   7 16

                            #autofilet_log = open('autofilet.log', 'a', 0)
                            autofilet_log = open('autofilet.log', 'a', 0)
                            autofilet_log.write('ERROR: %s, ext = %i xord = %i yord = %i was terminated because autofilet had not finished after 3 minutes, no output written. \n' %(ifile, extn, xord, yord))
                            
                            if (xord == 15) and (yord == 4):
                                xord = 17
                                yord = 4
                            elif (xord == 6) and (yord == 6):
                                xord = 7
                                yord = 7
                            elif (xord == 15) and (yord == 29):
                                xord = 15
                                yord = 29
                            elif (xord == 6) and (yord == 15):
                                xord = 6
                                yord = 17

                            elif (xord == 17) and (yord == 4):
                                xord = 17
                                yord = 5
                            elif (xord == 7) and (yord == 7):
                                xord = 6
                                yord = 8
                            elif (xord == 15) and (yord == 31):
                                xord = 17
                                yord = 29
                            elif (xord == 6) and (yord == 17):
                                xord = 7
                                yord = 15

                            elif (xord == 17) and (yord == 5):
                                xord = 18
                                yord = 5
                            elif (xord == 6) and (yord == 8):
                                xord = 7
                                yord = 8
                            elif (xord == 17) and (yord == 29):
                                xord = 16
                                yord = 30
                            elif (xord == 7) and (yord == 15):
                                xord = 7
                                yord = 16
                            elif (xord == 18 and yord == 5) or (xord == 7 and yord == 8) or (xord == 16 and yord == 30) or (xord == 7 and yord == 16):
                                print xord, yord
                                sys.exit('ERROR: %s, ext %i hung during calibration for 4 different xord and yord. Please remove from calibration or manually recalibrate' %(ifile, extn))                    
                            
        
                            autofilet_log.write('\tTrying again with xord = %i, yord = %i. \n' %(xord, yord))
                            autofilet_log.flush()
                            autofilet_log.close()
                            success = False
                        else :
                            success = True
                    autofilet_call.wait()

                
    #autofilet_log.close()
    inventory.close()
 
def build_inventory():
    '''
    The function parses the autofilet.log file to find where a file is not "VERIFIED"
    or where an ERROR occurred
    Inputs: none
    Outputs:
        rootname: the names of files
        extn: an array of science extension numbers corresponding to rootname
        warning: list of booleans as to whether an error was issued
    '''
    rootname = glob.glob('PSUB/*')
    for i, ifile in enumerate(rootname):
        rootname[i] = ifile.strip('PSUB/').split('_rwc')[0]
    rootname = np.array(rootname)
    ofile = open('autofilet.log', 'r')
    all_lines = ofile.readlines()
    extn = np.arange(len(rootname))
    warning = np.array([bool(warn) for warn in np.copy(rootname)])
    new_image_flag = True
    indx = 0
    new_run_indices = np.where(np.array(all_lines) == 'NEW RUN STARTS HERE\n')
    if len(new_run_indices[0]) < 1:
        new_run_indx = 0
    else:
        new_run_indx = new_run_indices[0][-1]
    log_dict = {}
    termination_flag = False
    for iline in all_lines[new_run_indx:]:
        if 'image=' in iline:
            new_image_flag = True
            tmp_rootname = iline.split('=')[1][5:-11].strip('_')  #get just the rootname
            tmp_ext = iline.split('=')[1].split('[')[1]
            tmp_warning = False #initialize to False, change if a warning is encountered
            termination_flag = False
        elif ('VERIFIED' in iline) & (new_image_flag == True):
            tmp_warning = True
            new_image_flag = False
        elif (('elapsed' in iline) & (not termination_flag)): #processing finished normally
            if tmp_rootname not in log_dict.keys():
                log_dict[tmp_rootname] = {'ext':[], 'warning':[]}
            log_dict[tmp_rootname]['ext'].append(tmp_ext)
            log_dict[tmp_rootname]['warning'].append(tmp_warning)
            warning[indx] = tmp_warning
            extn[indx] = tmp_ext
            indx += 1

        elif ('ERROR' in iline):
            print 'ERROR found in log file:'
            print '\t', iline
            termination_flag = True
    #warning = np.array([bool(warn) for warn in warning])
    #extn = np.array(extn)
    #rootname = np.array(rootname)
    if indx != len(rootname):
        print '!!!!!!!!!!!!!!!!\nWARNING, there are %i files in PSUB, but only %i are being rebuilt\n!!!!!!!!!!!!!!!!' %(len(rootname), indx)

    #pdb.set_trace()
    #return rootname, extn, warning
    return log_dict
            
       
def rebuild_files(interactive = False, failed_files = []):

    '''

    This code is used to rebuild the raw files from the PSUB files. 

    There are two cases of files which should not be automatically added to your final file
        1. Autofilet.pro was terminated before writing an output file
        2. Autofilet.pro created an output file, but was unsure it corrected the herringbone pattern correctly
    This code will never consider adding case 1.
    Case 2: If the code is being run interactively, then plots will be shown and the user will be asked to decide
            whether to keep the file. If the code is not being run interactively, any files the user wishes to skip should
            be specified in failed_files

    Inputs:
        interactive: run in interactive mode, this will cause data with warnings in the autofilet.log
                     file to be plotted and verified interactively. If interactive is set to false, it
                     will assume all data are ok. Default: False
        failed_files: if a list of failed files in PSUB is known, it may be passed in to be used if interactive = False
                Default = []
                I haven't tested this option
    Outputs:
        Creates herringbone corrected raw files in a CLEAN directory
        A list of failed files
        
    Calls to:
        build_inventory
        create_pdfs
        

    Known Issues:
        This code may cut off the last n in the comment for HPATPARS
        This code only rebuilds datasets which are newly processed as part of this run.
    '''
    #warn users that this is not being run interactively
    if not interactive:
        print 'WARNING: this script assumes that you have verified the correction in all images not listed in failed_files'

    #get information for autofilet.log
    #rootname, extn, warning = build_inventory()
    log_dict = build_inventory()
    #create a list of unique dataset names
    #uniq_rootname = list(set(rootname))
    #create output file directories
    if not os.path.exists('FAILED'):
        os.mkdir('FAILED')
    if not os.path.exists('CLEAN'):
        os.mkdir('CLEAN')
    #Open file to record list of files which will be determined to have failed fits
    failed_list = open('FAILED/failed_list.txt', 'w')
    #Loop over each dataset name (this is of the raw files)
    #for indiv_file in uniq_rootname:
    for indiv_file in log_dict.keys():
        print 'rebuilding ', indiv_file
        #Find all entries which correspond to that dataset (multiple imsets)
        #indx = np.where(rootname == indiv_file)[0]
        #sort_indx = np.argsort(extn[indx])
        #interactively check the pattern removal of all files which require verification from autofilet.log
        #If verified, then modify warning to keep them
        #if not verified, then put them on the failed list and copy them to the FAILED directory
        for i, warn, ext in zip(range(len(log_dict[indiv_file]['ext'])), log_dict[indiv_file]['warning'], log_dict[indiv_file]['ext']):
            ext = int(ext)
            print '\t extn ', ext
            #pdb.set_trace()
            if warn == True:  #a warning exists for this file
                if interactive == True:
                    #create plots
                    not_reviewed = create_pdfs(save = False, interactive = interactive, glob_str = 'PSUB/%s_rwc_%i.fits' %(indiv_file, ext))
                    if not_reviewed == False:
                        keep_file = raw_input('Is the herring bone pattern removed from this file? (y, n) ')
                        if keep_file == 'y':
                            log_dict[indiv_file]['warning'][i] = False  #reset warning if pattern removal is verified
                        else: #copy file to failed and add to failed list
                            shutil.copyfile('PSUB/%s_rwc_%i.fits' %(indiv_file, ext), 'FAILED/%s_rwc_%i.fits' %(indiv_file, ext))
                            failed_list.write('%s_rwc_%i.fits\n'  %(indiv_file, ext))

                    else:
                        log_dict[indiv_file]['warning'][i] = False
                    
                else: #if interactive is False check input list of failed files, and if it is not in the list, then reset then remove the warning

                    if ('%s_rwc_%i.fits' %(indiv_file, ext) not in failed_files) &  ('%s_wav_rwc_%i.fits' %(indiv_file, ext) not in failed_files):
                        log_dict[indiv_file]['warning'][i] = False
                    else: #copy file to failed and add to failed list
                        shutil.copyfile('PSUB/%s_rwc_%i.fits' %(indiv_file, ext), 'FAILED/%s_rwc_%i.fits' %(indiv_file, ext))
                        failed_list.write('%s_rwc_%i.fits\n'  %(indiv_file, ext))

        #--------------------------
        #Copy the original file to CLEAN and update the extensions with approved patterns
        #--------------------------
        #Copy the file and open in update mode
        
        try: #if the file is not a wavecal 
            shutil.copyfile('ORIG/%s_raw.fits' %(indiv_file), 'CLEAN/%s_raw.fits' %(indiv_file))
            ofile = pyfits.open('CLEAN/%s_raw.fits' %(indiv_file), mode = 'update') 
            wave = False
        except: #if the file is a wavecal
            pdb.set_trace()
            shutil.copyfile('ORIG/%s_wav.fits' %(indiv_file), 'CLEAN/%s_wav.fits' %(indiv_file))
            ofile = pyfits.open('CLEAN/%s_wav.fits' %(indiv_file), mode = 'update')
            wave = True
            
        #Loop over all imsets in each file
        ofile[0].header.set('HPATCORR', 'COMPLETE', 'remove herring-bone pattern noise', before = 'DQICORR')
        
        #for i, warn, ext in zip(range(len(indx)), warning[indx], extn[indx]):
        for i, warn, ext in zip(range(len(log_dict[indiv_file]['ext'])), log_dict[indiv_file]['warning'], log_dict[indiv_file]['ext']):
            #Add new header group to each extension
            ext = int(ext)
            ofile[ext].header.set('   ', '/HERRINGBONE PATTERN NOISE REMOVAL INFORMATION', after = 'MEANBLEV')
            ofile.flush()
            #If no warning, then update the data extension
            data = pyfits.getdata('PSUB/%s_rwc_%i.fits' %(indiv_file, ext), 0)
            psub_hist = pyfits.getheader('PSUB/%s_rwc_%i.fits' %(indiv_file, ext), 0).get_history()

            if warn == False:
                ofile[ext].data = data
                #FITS files can only have signed 16 bit integers. Signed integers go from -32768 to +32768.
                #raw files are contain unsigned integers (so we don't need the negative numbers) and we want
                #to be able to represent as high a count rate as possible. Unsigned 16 bit integers can go from
                # 0 to 65536. In order to get the full range of possible unsigned integers, OPUS subtracts 32768
                #from the raw image data and sets BZERO. The files written by AUTOFILET contain BITPIX = -32 (floating
                #point numbers) and therefore do not need BZERO set. This next section removed BZERO from the rebuilt 
                #files to avoid applying the shift to the data twice.
                try:
                    ofile[ext].header.remove('BZERO')
                    ofile[ext].header.remove('BSCALE')
                except ValueError:
                    print 'WARNING: No header keyword BZERO in %s[%i]' %(indiv_file, ext)
                

            #Add new header keywords to record herringbone processing parameters 
                           
            ofile[ext].header.set('HPATVER', psub_hist[2].strip(), 'program version')#, after = '   ')
            ofile[ext].header.set('HPATDATE', ' '.join(psub_hist[1].split()[3:]), 'date program was run', after = 'HPATVER')
            ofile[ext].header.set('HPATPARS', psub_hist[6].strip(), 'detected pattern', after = 'HPATDATE')
            if not warn:
                ofile[ext].header.set('HPATPASS', True, 'success (T), failed (F; frame is original*1.0)', after = 'HPATPARS') 
            else:
                ofile[ext].header.set('HPATPASS', False, 'success (T), failed (F; frame is original*1.0)', after = 'HPATPARS') 

        ofile.flush()
        ofile.close()
        
              

#To do:
#Write code to extract the names of files which fail execution - I'm not sure if I've done this
#Check the IDL termination


if __name__ == "__main__":
    '''
    This code replaces the directions in AAS_README in the stis2 folder from Rolf Jansen.
    It removed the herringbone pattern (fixed pattern readnoise) introduced with the move
    to side 2 electronics on STIS. This code is run in parts. The first part removed the 
    herringbone pattern from each STIS imset in the ORIG folder. The next part rebuilds the
    raw files with the herringbone corrected science data.

    Requirements:
        Folder called ORIG with original raw files
        Folder called PSUB for the corrected imsets
        Folder called TABLES
        herringbone_correction_confirmation_plots - a program for plotting if interactive = True
        autofilet.pro - IDL program written by Rolf Jansen, this can be run as a stand alone program
        rfunct.pro : called by autofilet.pro
        splflt.pro : called by autofilet.pro
        splinefit.pro : called by autofilet.pro

    Command Line Arguments:
        --dir : The directory where ORIG and PSUB exist if not the current directory. Default: current directory
        --ioff : Turn off interactive mode. This simply removes the confirmation of the removal
                of the herringbone pattern if autofilet is not sure if the patter was correctly
                removed. Default: interactive = True

    '''


    #Set command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', dest = 'directory', help = 'Enter the location of the ORIG and PSUB directories', default = '')
    parser.add_argument('--ioff', dest = 'interactive', action = 'store_false', help = 'Turn off interactive mode', default = True)
    directory = parser.parse_args().directory
    interactive = parser.parse_args().interactive
    if len(directory) < 2:
        directory = os.getcwd()
    os.chdir(directory)
    if not os.path.exists('TABLES'):
        os.mkdir('TABLES')
    if not os.path.exists('PSUB'):
        os.mkdir('PSUB')
    autofilet()
    rebuild_files(interactive = interactive)
