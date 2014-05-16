batch_run_autofilet was written to avoid running the successive shell scripts detailed in 
AA_README.txt found in the stis2 folder on http://stis2.sese.asu.edu/. This script avoids
dependencies on fortran and parsing error statements (which may vary from computer to computer).

It is a wrapper which allows you to run batched of data through autofilet. 

Original (uncorrected) raw files should be in the folder ORIG
The Herringbone subtracted single extension images will be put in the PSUB folder
The herringbone subtracted equivalent files to those in ORIG will be put in CLEAN
Files in which the herringbone correction was not removed will be put in the FAILED directory
Information on the pattern removed is in the TABLES directory

This code has not been tested as extensively as Rolf's code (the underlying autofilet.pro 
has not been modified). It was tested on imaging, spectroscopic, wavecal, flat field, bias
file, and dark file. Some of these files had multiple imsets.

I never encountered IDL hanging, there is code to stop the process, but it has never been tested.

All required directories and will be in the tar file. The code can be run by typing:
python batch_run_autofilet.py

The command line options are:
        --dir : The directory where ORIG and PSUB exist if not the current directory. Default: current directory
        --ioff : Turn off interactive mode. This simply removes the confirmation of the removal
                of the herringbone pattern if autofilet is not sure if the patter was correctly
                removed. Default: interactive = True

You must have IDL and python installed (including the packages numpy and pyfits).

Accuracy of headers and parameters has been verified. Due to the introduction of random
noise into the final product, data values were checked to make sure they were within 0.15
(approximate level of random noise added) of the results produced by running the shell script.
Files were verified by comparing the output of the shell commands with the output of this 
script using pyfits.FITSDiff with a tolerance of 0.15

The keywords BZERO and BSCALE have been removed from the output file since the data are
floating point.