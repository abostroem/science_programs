#! /usr/bin/env python

import sys
import os.path
import string
import math
import copy
from time import sleep
from Tkinter import *
import tkSimpleDialog
import tkFileDialog
import tkMessageBox

import numpy as N
try:
    from stsci import ndimage
except ImportError:
    import ndimage
try:
    from stsci.convolve import boxcar
except ImportError:
    from convolve import boxcar
import pyfits
from matplotlib import pylab
from matplotlib import cm
from matplotlib import ticker
from PIL import Image as PILImage

import celwcs

__version__ = "2012 November 6"

# reference image default min and max values
vmin = None
vmax = None

#opo image default max and min
vomin = None
vomax = None

# figure numbers for reference and public release images
REF_FIGURE = 1
OPO_FIGURE = 2

# indicates where subsampling was done for block averaging
LOCN_CORNER = 1
LOCN_CENTER = 2

# bin to approximately this many pixels in each axis
BINNED_SIZE = 2048.

RADIANStoDEGREES = 180. / math.pi

class ComputeWCS (object):

    def __init__ (self, root):

        self.root = root
        self.ref = Image (REF_FIGURE)   # Image object for reference image
        self.opo = Image (OPO_FIGURE)   # Image object for public release image

        self.angle = None               # Tkinter variable

        self.xdata = None               # mouse cursor X coordinate
        self.ydata = None               # mouse cursor Y coordinate

        self.total_points = 0           # total number of points marked
        self.npoints = 0                # excluding deleted points
        self.refXY = {}                 # index: x,y in reference image
        self.opoXY = {}                 # index: x,y in public release image
        self.residuals = {}             # index: x,y residuals of fit in ref

        # additional WCS info for the public release image
        self.opo_orientation = None             # degrees
        self.opo_scale = None                   # arcseconds per pixel
        # This is a measure of the skew of the public release image,
        # based on ratios of the CD matrix elements.
        self.opo_skew = None
        # The user may specify an alternate WCS ('A', 'B', etc.).
        self.altWCS_index = 0
        self.altWCS = None

        # If this flag is set to True, that means we need to write the
        # WCS info to stdout when the user quits.  It will be False if
        # the plate solution has not been computed yet or if the WCS info
        # has already been written to an output file.
        self.print_wcs = False                  # initial value

        # Plate solution
        self.x_coeff = None
        self.y_coeff = None
        self.x_inverse = None
        self.y_inverse = None

        self.frame = Frame (root)
        self.frame.grid()

        # Set interactive mode on, in case the user has 'interactive : False'
        # in the ~/.matplotlib/matplotlibrc file.
        pylab.ion()

        # Get the name of the reference image, and display it.
        dummy0 = Button (self.frame, text="display reference image",
                         command=self.refImage)
        dummy0.grid (row=0, column=0, sticky=W)

        # Get the name of the public release image, and display it.
        dummy0 = Button (self.frame, text="display public-release image",
                         command=self.opoImage)
        dummy0.grid (row=1, column=0, sticky=W)

        # include a dummy label to separate sections of the frame
        #dummy0 = Label (self.frame, text="")
        #dummy0.grid (row=2, column=0, sticky=W)

        # Option to set the min and max values and redisplay reference image.
        dummy0 = Button (self.frame, text="set reference image limits",
                         command=self.refLimits)
        dummy0.grid (row=2, column=0, sticky=W)

        # Option to set the min and max values and redisplay opo image.
        dummy0 = Button (self.frame, text="set public-release image limits",
                         command=self.opoLimits)
        dummy0.grid (row=3, column=0, sticky=W)


        # Option to rotate the reference image in 90-degree increments.
        rotframe = Frame (self.frame)
        self.angle = IntVar()
        self.angle.set (0)
        Radiobutton (rotframe, text="original orientation",
                     value=0, variable=self.angle,
                command=self.rotateRefImage).grid (row=0, column=0, sticky=W)
        Radiobutton (rotframe, text="90 degrees",
                     value=90, variable=self.angle,
                command=self.rotateRefImage).grid (row=1, column=0, sticky=W)
        Radiobutton (rotframe, text="180 degrees",
                     value=180, variable=self.angle,
                command=self.rotateRefImage).grid (row=2, column=0, sticky=W)
        Radiobutton (rotframe, text="270 degrees",
                     value=270, variable=self.angle,
                command=self.rotateRefImage).grid (row=3, column=0, sticky=W)
        rotframe.grid (row=4, column=0)

        # The QUIT button.
        dummy0 = Button (self.frame, text="QUIT", command=self.exitQuit)
        dummy0.grid (row=5, column=0, sticky=W)

        # Status line.
        self.status_label = Label (self.frame,
                text="Note:  display reference image", fg="red")
        self.status_label.grid (row=6, column=0, columnspan=2, sticky=W)

        # Add points.
        dummy1 = Button (self.frame, text="add points", command=self.addPoints)
        dummy1.grid (row=0, column=1, sticky=W)

        dummy1 = Button (self.frame, text="stop adding points",
                         command=self.stopAddingPoints)
        dummy1.grid (row=1, column=1, sticky=W)

        # Add points using the fit.
        dummy1 = Button (self.frame, text="add points using fit",
                         command=self.addPointsUsingFit)
        dummy1.grid (row=2, column=1, sticky=W)

        # Delete a point.
        dummy1 = Button (self.frame, text="delete a point",
                         command=self.delPoint)
        dummy1.grid (row=3, column=1, sticky=W)

        # Public-release image may be negative.
        dummy1 = Button (self.frame, text="invert public-release image",
                         command=self.opoInvert)
        dummy1.grid (row=4, column=1, sticky=W)

        # Save coordinate info.
        dummy1 = Button (self.frame, text="write wcs info to file",
                         command=self.writeWCS)
        dummy1.grid (row=5, column=1, sticky=W)

        # print info about using wcs
        dummy2 = Button (self.frame, text="help", command=self.helpInfo)
        dummy2.grid (row=0, column=2, sticky=W)

    def exitQuit (self):
        if self.print_wcs:
            # Coordinate information has not been saved to a file.
            # Do you want to save it?
            save_them = tkMessageBox.askyesno ("parameters not saved",
                    "Coordinate parameters have not been saved to a file.\n" \
                    "Do you want to save them?", default="yes")
            if save_them:
                self.writeWCS()
                self.print_wcs = False
            if self.print_wcs:
                self.writeInfo (sys.stdout)
        self.frame.quit()

    def helpInfo (self):
        """Print info about using wcs.py."""

        print (
"""button "display reference image":
The reference image and public-release image must be specified before other
operations can be done.
Use this button to specify the name of the FITS-format reference image.
""")

        print (
"""button "display public-release image":
Use this button to specify the name of the public-release image (the one
for which coordinate information is to be computed).  The file format may be
FITS, TIFF, or JPG.
""")

        print (
"""button "set reference image limits":
By default the image display includes the full range of data values in the
input.  "set reference image limits" can be used to set the minimum and
maximum values for the reference image display.
""")

        print (
"""buttons "original orientation", "90 degrees", "180 degrees", "270 degrees":
If the reference image and public-release image are rotated significantly
with respect to each other, it can be difficult to match stars in the two
images.  This set of "radio" buttons lets you rotate the reference image by
a multiple of 90 degrees in order to approximately match the orientations
of the two images.
""")

        print (
"""button "QUIT":
Exit the application, closing all windows.
If the coordinate information has not been written to a file (see "write wcs
info to file") and you choose not to do so at this time, the info will be
written to the standard output.
""")

        print (
"""button "add points":
Click on this button to enter the mode for marking matching stars in the
two images.
Find a star (or other sharp feature) that is visible in both images.
Click on the star in the reference image, and then click on the star in
the public-release image.  Repeat.  In each image, the brightest pixel
near the cursor (with quadratic weighting centered on the cursor position)
will be taken as the location of the star.  Exit this mode with "stop adding
points" or by right-clicking in the reference image.
After four or more pairs of matching stars have been marked, a fit will
be computed, and the residuals of the fit will be printed to the standard
output.  Large "skew factors" or large residuals imply that one or more
points are mismatched (see "delete a point").
""")

        print (
"""button "stop adding points":
Click on this button to exit the "add points" mode.  However, due to a
bug in the code, you must also click once more in the reference image
before the program will actually leave the "add points" mode.
An alternative to using this button is to click in the reference image
with the right mouse button, which will exit the "add points" mode
immediately.
""")

        print (
"""button "add points using fit":
After marking at least three pairs of matching stars, you can use this
button to enter an "add points" mode in which you only mark stars in the
reference image.  The current fit between the two images will then be
used to compute the expected position of the matching star in the
public-release image, and the brightest pixel near that location (with
quadratic weighting, as usual) will be taken as the location of the
matching star in the public-release image.
You can exit this mode by using "stop adding points" or by right-clicking
with the mouse in the reference image.
""")

        print (
"""button "delete a point":
If it looks as if a pair of stars has been marked incorrectly (i.e. the
points in the two images aren't really the same star), as implied by the
residuals of the fit, that point can be deleted.  First exit "add points"
mode, then click "delete a point", then click on the point in the reference
image.  The images will be redisplayed.  Note that the numbers attached to
the other points will remain unchanged.
""")

        print (
"""button "invert public-release image":
Sometimes a public-release image may be a negative image (though it may be
displayed as a positive image due to the look-up table).  In that case when
you try to mark points ("add points") in that image, the plus sign showing
the location will not be centered on the star; this is because the program
looks for the brightest pixel, but the location would be a minimum in a
negative image.  Click on this button to invert the image in memory.
""")

        print (
"""button "write wcs info to file":
Click on this button to specify the name of a (new) text file to which the
WCS information will be written.  If this is not done, you will also have
the opportunity to specify a file name after clicking on "QUIT" (and if you
still do not specify a file, the information will be written to the screen).
""")

        if self.ref.data is None:
            self.showStatus ("Note:  display reference image")
        elif self.opo.data is None:
            self.showStatus ("Note:  display public-release image")
        else:
            self.showStatus ("")

    def showStatus (self, message):
        """Display message on the status line."""
        self.status_label.config (text=message)
        self.status_label.update_idletasks()

    def refImage (self):
        """Get the name of the reference file, and display the image."""
        self.ref.name = tkFileDialog.askopenfilename (
                 title="reference image", parent=self.frame)
        if not self.ref.name:
            return
        self.showStatus ("reading reference image ...")

        # Reset these, in case another image was loaded earlier.
        self.total_points = 0
        self.npoints = 0
        self.refXY = {}
        self.opoXY = {}
        self.residuals = {}
        self.angle.set (0)

        (tempdata, self.ref.header) = pyfits.getdata (self.ref.name,
                header=True)
        shape = tempdata.shape
        n = math.sqrt (shape[0] * shape[1])
        block_size = int (round (float (n) / BINNED_SIZE))
        block_size = max (1, block_size)
        if block_size > 1:
            self.ref.bin_locn = LOCN_CORNER
            boxcar (tempdata, (block_size,block_size), output=tempdata)
            self.ref.data = tempdata[0:-1:block_size,0:-1:block_size]
            del tempdata
            print "info:  binning for reference image is", block_size
        else:
            self.ref.data = tempdata
        self.ref.current_data = self.ref.data

        # Set binning factor, compute mapping to binned pixels.
        self.ref.setBinning (block_size)

        if self.ref.Vmin is None:
            self.ref.Vmin = N.minimum.reduce (self.ref.data.flat)
        if self.ref.Vmax is None:
            self.ref.Vmax = N.maximum.reduce (self.ref.data.flat)

        self.displayRefImage()
        (self.ref.xlow, self.ref.xhigh) = pylab.xlim()
        (self.ref.ylow, self.ref.yhigh) = pylab.ylim()

        shape = self.ref.data.shape
        self.ref.marklen = self.ref.length_fraction * max (shape[0], shape[1])
        self.ref.setSearchParameters()

        if self.opo.data is None:
            self.showStatus ("Note:  display public-release image")
        else:
            self.showStatus ("")

        # Read the WCS info from the reference image header.
        self.getRefWCS()

    def getRefWCS (self):
        """Get the WCS info for the reference image."""

        hdr = self.ref.header
        shape = self.ref.data.shape
        # shape of the binned data
        self.ref.bnaxis[0] = shape[1]           # note:  swapping axis numbers
        self.ref.bnaxis[1] = shape[0]
        # shape of the original data
        self.ref.naxis[0] = hdr["naxis1"]
        self.ref.naxis[1] = hdr["naxis2"]
        no_wcs = False                          # initial value
        alt = [""]
        if hdr.has_key ("CRPIX1"):              # likely to be present
            if hdr.has_key ("WCSNAME"):
                wcsname = [hdr["WCSNAME"]]
            else:
                wcsname = ["primary WCS"]
            for a in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
                      'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R',
                      'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']:
                if hdr.has_key ("CRPIX1" + a):
                    alt.append (a)
                    key = "WCSNAME" + a
                    if hdr.has_key (key):
                        wcsname.append (hdr[key])
                    else:
                        wcsname.append ("")
            if len (alt) > 1:
                # ask which WCS to use
                self.inputAlt (alt, wcsname)
        try:
            # read WCS info from header
            self.altWCS = alt[self.altWCS_index]
            self.ref.wcs = celwcs.celwcs (hdr, alternate=self.altWCS)
            (self.ref.crval, self.ref.crpix, self.ref.cd) = \
                    self.ref.wcs.getRADec_WCS()
            self.ref.current_crpix = self.ref.crpix.copy()
            self.ref.current_cd = self.ref.cd.copy()
        except ValueError:
            no_wcs = True
        if no_wcs:
            if self.altWCS:
                tkMessageBox.showerror ("no WCS info",
                    "Alternate WCS '%s' not found in reference image %s" % \
                    (self.altWCS, self.ref.name))
            else:
                tkMessageBox.showerror ("no WCS info",
                    "WCS info not found in reference image %s" % self.ref.name)
            pylab.close (REF_FIGURE)
            self.ref.data = None
            self.showStatus ("Note:  redisplay reference image")
        elif len (alt) > 1:
            if self.altWCS:
                print "info:  using alternate WCS '%s'" % self.altWCS
            else:
                print "info:  using primary WCS"

    def inputAlt (self, alt, wcsname):
        """Let the user specify which WCS to read.

        @param alt: list of WCS identifiers, "" for the primary or a single
            letter for an alternate WCS
        @type alt: list

        @param wcsname: list of the name (or a null string) for each WCS
        @type wcsname: list of strings
        """

        self.alt_top = Toplevel (self.frame)
        altframe = Frame (self.alt_top)
        altframe.grid()
        txt = Label (altframe, text="specify which WCS to use")
        txt.grid (row=0, column=0, sticky=W)
        scrollbar = Scrollbar (altframe, orient=VERTICAL)
        scrollbar.grid (row=1, column=1, sticky=NS)
        self.listbox = Listbox (altframe, height=15,
                                selectmode=SINGLE,
                                yscrollcommand=scrollbar.set)
        self.listbox.grid (row=1, column=0, sticky=W)
        scrollbar.config (command=self.listbox.yview)
        for i in range (len (alt)):
            self.listbox.insert (END, "%1s  %s" % (alt[i], wcsname[i]))
        self.listbox.bind ("<Button1-ButtonRelease>", self.setAlt)
        self.root.wait_window (self.alt_top)

        self.showStatus ("")

    def setAlt (self, event):

        value = self.listbox.curselection()
        if len (value) < 1:
            self.altWCS_index = ""
        self.altWCS_index = int (value[0])
        self.alt_top.destroy()

    def opoImage (self):
        """Get the name of the public release file, and display the image."""
        self.opo.name = tkFileDialog.askopenfilename (
                 title="public-release image", parent=self.frame)
        if not self.opo.name:
            return
        self.showStatus ("reading public-release image ...")

        # Reset these, in case another image was loaded earlier.
        self.total_points = 0
        self.npoints = 0
        self.refXY = {}
        self.opoXY = {}
        self.residuals = {}

        if self.opo.name.endswith (".fits") or \
           self.opo.name.endswith (".fit"):
            tempdata = pyfits.getdata (self.opo.name)
            shape = tempdata.shape
            self.opo.crpix[0] = shape[1] / 2.
            self.opo.crpix[1] = shape[0] / 2.
            n = math.sqrt (shape[0] * shape[1])
            block_size = int (round (float (n) / BINNED_SIZE))
            block_size = max (1, block_size)
            if block_size > 1:
                self.opo.bin_locn = LOCN_CORNER
                boxcar (tempdata, (block_size,block_size), output=tempdata)
                self.opo.data = tempdata[0:-1:block_size,0:-1:block_size]
                del tempdata
                print "info:  binning for public-release image is", block_size
            else:
                self.opo.data = tempdata
        else:
            im = PILImage.open (self.opo.name)
            size = im.size
            self.opo.crpix[0] = size[0] / 2.
            self.opo.crpix[1] = size[1] / 2.
            if len (size) != 2:
                raise RuntimeError, "The shape of %s is peculiar" % \
                        self.opo.name
            if im.mode != "RGB":
                tkMessageBox.showwarning ("unknown mode",
                    "the mode %s of %s might be a problem" % \
                            (im.mode, self.opo.name))
            n = math.sqrt (size[0] * size[1])
            block_size = int (round (float (n) / BINNED_SIZE))
            block_size = max (1, block_size)
            if block_size > 1:
                self.opo.bin_locn = LOCN_CENTER
                # make the block size an odd number
                if (block_size // 2 * 2) == block_size:
                    block_size += 1
                print "info:  size of public-release image is", size
                print "info:  binning for public-release image is", block_size
                oposize = [size[0] // block_size, size[1] // block_size]
                # round up the X axis length (but not the Y axis)
                if oposize[0] * block_size < size[0]:
                    oposize[0] += 1
                if oposize[0] < 1 or oposize[1] < 1:
                    tkMessageBox.showwarning ("too small",
                        "can't handle public-release image %s" % self.opo.name)
                    print "ERROR:  original size is", size
                    print "ERROR:  binned size is", oposize
                    if oposize[0] < 1:
                        oposize[0] = 1
                    if oposize[1] < 1:
                        oposize[1] = 1
                self.opo.data = N.zeros ((oposize[1], oposize[0], 3),
                                         dtype=N.float32)
                ylow = 0
                yhigh = block_size
                half_block = block_size // 2
                for n in range (oposize[1]):
                    y = oposize[1] - 1 - n
                    if yhigh > size[1]:
                        yhigh = size[1]
                        if yhigh - ylow < block_size:
                            break
                    box = (0, ylow, size[0], yhigh)
                    region = im.crop (box)
                    regionstr = region.tostring()
                    tempdata = N.fromstring (regionstr, dtype=N.uint8)
                    tempdata = N.reshape (tempdata,
                                          newshape=(yhigh-ylow, size[0], 3))
                    nregion = tempdata.astype (N.float32) / 255.
                    # block average each color separately;
                    temp0 = boxcar (nregion[:,:,0], (block_size,block_size))
                    temp1 = boxcar (nregion[:,:,1], (block_size,block_size))
                    temp2 = boxcar (nregion[:,:,2], (block_size,block_size))
                    # fill opo.data from the top down
                    if n == 0:
                        # this is a 1-D subarray
                        nelem = temp0[0,half_block:-1:block_size].shape[0]
                    # sample at the center of each block
                    self.opo.data[y,0:nelem,0] = \
                                temp0[half_block,half_block:-1:block_size]
                    self.opo.data[y,0:nelem,1] = \
                                temp1[half_block,half_block:-1:block_size]
                    self.opo.data[y,0:nelem,2] = \
                                temp2[half_block,half_block:-1:block_size]
                    ylow += block_size
                    yhigh += block_size
                del tempdata
            else:
                print "info:  size of public-release image is", size
                imstr = im.tostring()
                del im
                nim = N.fromstring (imstr, dtype=N.uint8)
                nim = N.reshape (nim, newshape=(size[1], size[0], 3))
                del imstr
                self.opo.data = nim[::-1,:,:].astype (N.float32) / 255.
                del nim

        # set binning factor, compute mapping to binned pixels.
        self.opo.setBinning (block_size)
        if self.opo.Vmin is None:
            self.opo.Vmin = N.minimum.reduce (self.opo.data.flat)
        if self.opo.Vmax is None:
            self.opo.Vmax = N.maximum.reduce (self.opo.data.flat)
        self.displayOPOImage()
        (self.opo.xlow, self.opo.xhigh) = pylab.xlim()
        (self.opo.ylow, self.opo.yhigh) = pylab.ylim()

        shape = self.opo.data.shape
        self.opo.marklen = self.opo.length_fraction * max (shape[0], shape[1])
        self.opo.setSearchParameters()

        # Mark the center of the image.
        self.opo.drawPlus (x0=self.opo.crpix[0], y0=self.opo.crpix[1],
                           lenfactor=2, color="r", rescale=True)

        if self.ref.data is None:
            self.showStatus ("Note:  display reference image")
        else:
            self.showStatus ("")

    def refLimits (self):
        global vmin, vmax
        vmin = self.ref.Vmin
        vmax = self.ref.Vmax
        x = GetRefMinMax (self.frame, title="give min and max")
        if x.vmin is None or x.vmax is None:
            return
        self.ref.Vmin = x.vmin
        self.ref.Vmax = x.vmax
        if self.ref.data is None:
            self.refImage()
        else:
            self.displayRefImage()
            if self.npoints > 0:
                for i in range (1, self.total_points + 1):
                    if self.refXY.has_key (i):
                        (xr, yr) = self.refXY[i]
                        self.ref.drawPlus (xr, yr, index=i, color="c",
                                           draw_it=False)
                pylab.draw()
                #pylab.xlim (self.ref.xlow, self.ref.xhigh)
                #pylab.ylim (self.ref.ylow, self.ref.yhigh)


    def opoLimits (self):
        global vomin, vomax
        vomin = self.opo.Vmin
        vomax = self.opo.Vmax
        x = GetOpoMinMax (self.frame, title="give min and max")
        #x.vomin = 0
        #x.vomax = 10
        if x.vomin is None or x.vomax is None:
            return
        self.opo.Vmin = x.vomin
        self.opo.Vmax = x.vomax
        if self.opo.data is None:
            self.opoImage()
        else:
            self.displayOPOImage()
            if self.npoints > 0:
                for i in range (1, self.total_points + 1):
                    if self.opoXY.has_key (i):
                        (xr, yr) = self.opoXY[i]
                        self.opo.drawPlus (xr, yr, index=i, color="c",
                                           draw_it=False)
                pylab.draw()
                #pylab.xlim (self.opo.xlow, self.opo.xhigh)
                #pylab.ylim (self.opo.ylow, self.opo.yhigh)

    def rotateRefImage (self):
        if self.ref.data is None:
            self.refImage()
            if self.ref.data is None:
                tkMessageBox.showerror ("no reference image",
                "you must specify a reference image before you can rotate it")
                return
        self.ref.rotangle = int (self.angle.get())
        angle = float (self.ref.rotangle)
        self.showStatus ("rotating reference image ...")
        if self.ref.rotangle == 0:
            self.ref.current_data = self.ref.data
            self.ref.current_crpix = self.ref.crpix.copy()
            self.ref.current_cd = self.ref.cd.copy()
        else:
            self.ref.current_data = ndimage.rotate (self.ref.data, angle)
            if self.ref.rotangle == 90:
                self.ref.current_crpix[0] = self.ref.crpix[1]
                self.ref.current_crpix[1] = \
                        self.ref.bnaxis[0] * self.ref.block_size - \
                        self.ref.crpix[0] - 1.
                cosa = 0.
                sina = 1.
            elif self.ref.rotangle == 180:
                self.ref.current_crpix[0] = \
                        self.ref.bnaxis[0] * self.ref.block_size - \
                        self.ref.crpix[0] - 1.
                self.ref.current_crpix[1] = \
                        self.ref.bnaxis[1] * self.ref.block_size - \
                        self.ref.crpix[1] - 1.
                cosa = -1.
                sina = 0.
            elif self.ref.rotangle == 270:
                self.ref.current_crpix[0] = \
                        self.ref.bnaxis[1] * self.ref.block_size - \
                        self.ref.crpix[1] - 1.
                self.ref.current_crpix[1] = self.ref.crpix[0]
                cosa = 0.
                sina = -1.
            else:
                print "internal error, angle %d not supported" % \
                      self.ref.rotangle
            m_rot = N.matrix (((cosa, -sina), (sina, cosa)), dtype=N.float64)
            self.ref.current_cd = (N.matrix (self.ref.cd) * m_rot).A

        # Update the celwcs object with the current crpix and CD matrix.
        self.ref.wcs.setRADec_WCS (crpix=self.ref.current_crpix,
                                   cd=self.ref.current_cd)

        self.displayRefImage()
        if self.opo.data is None:
            self.showStatus ("Note:  display public-release image")
        else:
            self.showStatus ("")

        # xxx It would be better to modify these rather than clobbering them.
        if self.total_points > 0:
            self.total_points = 0
            self.npoints = 0
            self.refXY = {}
            self.opoXY = {}
            if self.opo.data is not None:
                self.displayOPOImage()
                self.opo.drawPlus (x0=self.opo.crpix[0], y0=self.opo.crpix[1],
                                   lenfactor=2, color="r", draw_it=False)
                pylab.xlim (self.opo.xlow, self.opo.xhigh)
                pylab.ylim (self.opo.ylow, self.opo.yhigh)

    def opoInvert (self):
        """Invert the public release image."""

        if self.opo.data is not None:
            self.opo.data = 1. - self.opo.data
            self.refreshImageDisplays()

    def addPoints (self):

        self.addPointsDone = False
        self.showStatus ("Mark points, reference first")

        pylab.figure (REF_FIGURE)
        (self.ref.xlow, self.ref.xhigh) = pylab.xlim()
        (self.ref.ylow, self.ref.yhigh) = pylab.ylim()
        pylab.figure (OPO_FIGURE)
        (self.opo.xlow, self.opo.xhigh) = pylab.xlim()
        (self.opo.ylow, self.opo.yhigh) = pylab.ylim()

        while not self.addPointsDone:

            self.read_mouse_xy (REF_FIGURE)
            if self.addPointsDone:
                break
            (x0_b, y0_b) = self.locatePoint (self.ref.current_data,
                               self.xdata, self.ydata,
                               search=self.ref.search, weight=self.ref.weight)
            if x0_b is None:
                tkMessageBox.showinfo ("try again",
                    "the cursor was outside the image\ntry again")
                self.showStatus ("add points:  mark in reference image")
                continue
            else:
                (x0, y0) = self.ref.fromBinned ((x0_b, y0_b))
                self.refXY[self.total_points+1] = (x0, y0)
                self.ref.drawPlus (x0, y0, index=self.total_points+1,
                                   color="c")
                self.showStatus ("add points:  mark in public-release image")
            pylab.xlim (self.ref.xlow, self.ref.xhigh)
            pylab.ylim (self.ref.ylow, self.ref.yhigh)

            self.read_mouse_xy (OPO_FIGURE)
            (x0_b, y0_b) = self.locatePoint (self.opo.data,
                               self.xdata, self.ydata,
                               search=self.opo.search, weight=self.opo.weight)
            while x0_b is None:
                tkMessageBox.showinfo ("try again",
                    "the cursor was outside the image\ntry again")
                self.showStatus ("add points:  mark in public-release image")
                self.read_mouse_xy (OPO_FIGURE)
                (x0_b, y0_b) = self.locatePoint (self.opo.data,
                               self.xdata, self.ydata,
                               search=self.opo.search, weight=self.opo.weight)
            (x0, y0) = self.opo.fromBinned ((x0_b, y0_b))
            self.opoXY[self.total_points+1] = (x0, y0)
            self.opo.drawPlus (x0, y0, index=self.total_points+1, color="c")
            self.showStatus ("add points:  mark in reference image")
            pylab.xlim (self.opo.xlow, self.opo.xhigh)
            pylab.ylim (self.opo.ylow, self.opo.yhigh)

            self.total_points += 1
            self.npoints += 1

            if self.npoints >= 3:
                # Do a plate solution, compute residuals, compute WCS.
                self.computeOpoWCS()
                self.showResiduals()

        self.showStatus ("")

    def stopAddingPoints (self):

        # xxx This isn't sufficient; I have to click one additional time in
        # the reference image window.  How can this be fixed?
        self.addPointsDone = True
        self.showStatus ("click once more in reference image (anywhere)")

    def addPointsUsingFit (self):

        if self.npoints < 3:
            tkMessageBox.showwarning ("not enough points",
                "you must mark at least three points before using this option")
            self.showStatus ("")
            return

        self.addPointsDone = False
        self.showStatus ("Mark points in reference image")

        pylab.figure (REF_FIGURE)
        (self.ref.xlow, self.ref.xhigh) = pylab.xlim()
        (self.ref.ylow, self.ref.yhigh) = pylab.ylim()
        pylab.figure (OPO_FIGURE)
        (self.opo.xlow, self.opo.xhigh) = pylab.xlim()
        (self.opo.ylow, self.opo.yhigh) = pylab.ylim()

        while not self.addPointsDone:

            self.read_mouse_xy (REF_FIGURE)
            if self.addPointsDone:
                break
            # _b in a variable name for a coordinate means that the
            # coordinate is in binned pixels
            (x0_b_ref, y0_b_ref) = self.locatePoint (self.ref.current_data,
                           self.xdata, self.ydata,
                           search=self.ref.search, weight=self.ref.weight)
            if x0_b_ref is None:
                tkMessageBox.showinfo ("outside image",
                    "the cursor was outside the image\ntry again")
                continue
            (x0_ref, y0_ref) = self.ref.fromBinned ((x0_b_ref, y0_b_ref))
            # this is the predicted position in the public-release image
            (x_fit, y_fit) = self.refToOpo (x0_ref, y0_ref)
            if x_fit is None:
                self.showStatus ("")
                return
            # predicted position in public-release image, but in binned pixels
            (x_b_fit, y_b_fit) = self.opo.toBinned ((x_fit, y_fit))
            (x0_b_opo, y0_b_opo) = self.locatePoint (self.opo.data,
                       x_b_fit, y_b_fit,
                       search=self.opo.search, weight=self.opo.weight)
            if x0_b_opo is None:
                tkMessageBox.showinfo ("outside image",
                    "the corresponding point is outside " \
                    "the public-release image\ntry again")
                continue
            (x0_opo, y0_opo) = self.opo.fromBinned ((x0_b_opo, y0_b_opo))
            self.refXY[self.total_points+1] = (x0_ref, y0_ref)
            self.opoXY[self.total_points+1] = (x0_opo, y0_opo)
            self.ref.drawPlus (x0_ref, y0_ref,
                               index=self.total_points+1, color="c")
            pylab.xlim (self.ref.xlow, self.ref.xhigh)
            pylab.ylim (self.ref.ylow, self.ref.yhigh)
            self.opo.drawPlus (x0_opo, y0_opo,
                               index=self.total_points+1, color="c")
            pylab.xlim (self.opo.xlow, self.opo.xhigh)
            pylab.ylim (self.opo.ylow, self.opo.yhigh)

            self.total_points += 1
            self.npoints += 1

            if self.npoints >= 3:
                # Do a plate solution, compute residuals, compute WCS.
                self.computeOpoWCS()
                self.showResiduals()

        self.showStatus ("")

    def delPoint (self):

        if self.npoints < 1:
            tkMessageBox.showinfo ("no points",
                    "There are no marked points to delete.")
            return

        pylab.figure (REF_FIGURE)
        (self.ref.xlow, self.ref.xhigh) = pylab.xlim()
        (self.ref.ylow, self.ref.yhigh) = pylab.ylim()
        pylab.figure (OPO_FIGURE)
        (self.opo.xlow, self.opo.xhigh) = pylab.xlim()
        (self.opo.ylow, self.opo.yhigh) = pylab.ylim()

        self.showStatus ("Mark point in reference image")

        # Note that we work with binned pixel coordinates here.
        self.read_mouse_xy (REF_FIGURE)
        nfound = 0
        mindist = None
        index = -1
        for i in range (1, self.total_points + 1):
            if self.refXY.has_key (i):
                (xr, yr) = self.ref.toBinned (self.refXY[i])
                distance = math.sqrt ((xr - self.xdata)**2 +
                                      (yr - self.ydata)**2)
                if distance <= self.ref.search:
                    nfound += 1
                    if mindist is None or distance < mindist:
                        mindist = distance
                        index = i
        if nfound > 0:
            if nfound > 1:
                tkMessageBox.showwarning ("ambiguous",
                        "deleting one of %d nearby points" % nfound)
            self.delThisIndex (index)
        else:
            tkMessageBox.showinfo ("not found",
                    "No nearby marked point found in reference image.")
            self.showStatus ("")

    def delThisIndex (self, index):
        """Delete the point with the specified index."""

        del (self.refXY[index])
        del (self.opoXY[index])
        if self.residuals.has_key (index):
            del (self.residuals[index])
        self.npoints -= 1
        self.showStatus ("point deleted (refreshing displays) ...")
        self.refreshImageDisplays()
        self.computeOpoWCS()                    # recompute with fewer points
        self.showResiduals()
        self.showStatus ("")

    def refreshImageDisplays (self):

        self.displayRefImage()
        self.displayOPOImage()
        self.opo.drawPlus (x0=self.opo.crpix[0], y0=self.opo.crpix[1],
                           lenfactor=2, color="r", draw_it=False)

        for i in range (1, self.total_points + 1):
            if self.refXY.has_key (i):
                (xr, yr) = self.refXY[i]
                self.ref.drawPlus (xr, yr, index=i, color="c", draw_it=False)
        pylab.draw()
        pylab.xlim (self.ref.xlow, self.ref.xhigh)
        pylab.ylim (self.ref.ylow, self.ref.yhigh)
        for i in range (1, self.total_points + 1):
            if self.opoXY.has_key (i):
                (xo, yo) = self.opoXY[i]
                self.opo.drawPlus (xo, yo, index=i, color="c", draw_it=False)
        pylab.draw()
        pylab.xlim (self.opo.xlow, self.opo.xhigh)
        pylab.ylim (self.opo.ylow, self.opo.yhigh)

    def showResiduals (self):

        if self.npoints < 3:
            return

        if self.residuals.has_key ("rms_x"):
            if self.npoints > 3:
                rms = math.sqrt (self.residuals["rms_x"]**2 +
                                 self.residuals["rms_y"]**2)
                print "RMS = %.3g; skew factor = %.3f" % (rms, self.opo_skew)
            else:
                print "skew factor = %.3f" % self.opo_skew
        if self.npoints > 3:
            print "residuals of fit:"
            for i in range (1, self.total_points + 1):
                if self.residuals.has_key (i):
                    (x_resid, y_resid) = self.residuals[i]
                    print "%d:  %5.2f" % \
                        (i, math.sqrt (x_resid**2 + y_resid**2))

    def locatePoint (self, image, x0, y0, search=10, weight=None):
        """Find the brightest point near x0, y0."""

        # Note that these are binned pixel coordinates.
        shape = image.shape
        nx = shape[1]
        ny = shape[0]
        x0 = int (round (x0))
        y0 = int (round (y0))
        x1 = x0 - search
        y1 = y0 - search
        x2 = x1 + 2*search+1
        y2 = y1 + 2*search+1

        if x1 < 0 or y1 < 0 or x2 >= nx or y2 >= ny:
            return (None, None)

        region = image[y1:y2,x1:x2].copy().astype (N.float32)
        # Use quadratic weighting centered on the user-specified
        # location to reduce the chance of finding a bright point
        # that is close to but offset from the intended target.
        region *= weight
        position = ndimage.maximum_position (region)
        x = position[1] + x1
        y = position[0] + y1

        return (x, y)

    def computeOpoWCS (self):

        if self.opo.wcs is None:
            # Assign initial values.
            self.opo.wcs = celwcs.celwcs ({"wcsdim": 2,
                    "crpix1": self.opo.crpix[0], "crpix2": self.opo.crpix[1],
                    "crval1": 0., "crval2": 0.,
                    "cdelt1": 1., "cdelt2": 1.,
                    "pc1_1": 1., "pc1_2": 0., "pc2_1": 0., "pc2_2": 1.,
                    "ctype1": "RA---TAN", "ctype2": "DEC--TAN"})

        if self.npoints < 3:
            self.print_wcs = False
            self.opo.crval[0] = None
            self.opo.crval[1] = None
            self.opo_orientation = None
            self.opo_scale = None
            return

        self.x_coeff, self.y_coeff = self.plateSolution()

        self.computeResiduals()

        # Find the point in the reference image that corresponds to the
        # center of the public release image, and use the WCS in the former
        # to convert to tangent-plane coordinates.
        (xr0, yr0) = self.opoToRef (self.opo.crpix[0], self.opo.crpix[1])
        (xi0, eta0) = self.xi_eta (xr0, yr0)

        crval = self.ref.wcs.frompixel ((xr0, yr0))
        self.opo.crval[0] = crval[self.ref.wcs.longitudeAxis()]
        self.opo.crval[1] = crval[self.ref.wcs.latitudeAxis()]

        # This is the rather arbitrary length of a vector in the +X or +Y
        # direction, to reduce numerical error when subtracting.
        v_length = 1000.

        # _x and _y indicate offsets in x and y respectively
        (xr_x, yr_x) = self.opoToRef (self.opo.crpix[0] + v_length,
                                      self.opo.crpix[1])
        (xi_x, eta_x) = self.xi_eta (xr_x, yr_x)

        (xr_y, yr_y) = self.opoToRef (self.opo.crpix[0],
                                      self.opo.crpix[1] + v_length)
        (xi_y, eta_y) = self.xi_eta (xr_y, yr_y)

        self.opo.cd[0,0] = (xi_x - xi0) / v_length      # d(xi) / dx
        self.opo.cd[0,1] = (xi_y - xi0) / v_length      # d(xi) / dy
        self.opo.cd[1,0] = (eta_x - eta0) / v_length    # d(eta) / dx
        self.opo.cd[1,1] = (eta_y - eta0) / v_length    # d(eta) / dy

        # arcseconds per pixel in the public release image
        self.opo_scale = math.sqrt ((xi_y - xi0)**2 + (eta_y - eta0)**2) * \
                3600. / v_length
        # orientation of Y axis, degrees eastward from north
        self.opo_orientation = math.atan2 (xi_y - xi0, eta_y - eta0) * \
                RADIANStoDEGREES

        # this is an attempt to measure the skew
        if abs (self.opo.cd[0,0]) * abs (self.opo.cd[1,1]) >= \
           abs (self.opo.cd[0,1]) * abs (self.opo.cd[1,0]):
            self.opo_skew = -self.opo.cd[1,1] / self.opo.cd[0,0]
        else:
            self.opo_skew =  self.opo.cd[0,1] / self.opo.cd[1,0]

        # Copy the WCS parameters into the celwcs object.
        self.opo.wcs.setRADec_WCS (self.opo.crval, self.opo.crpix, self.opo.cd)

        # Set this to true because we have computed new WCS info.
        self.print_wcs = True

    def plateSolution (self):
        """Get the mapping from public release image to reference image."""

        if self.npoints < 3:
            return None, None

        matrix = N.identity (3, dtype=N.float64)
        x_vector = N.zeros (3, dtype=N.float64)
        y_vector = N.zeros (3, dtype=N.float64)
        self.computeSums (matrix, x_vector, y_vector)

        x_coeff = N.linalg.solve (matrix, x_vector)
        y_coeff = N.linalg.solve (matrix, y_vector)

        return x_coeff, y_coeff

    def computeSums (self, matrix, x_vector, y_vector):
        """Compute sums for a linear equation.

        These arguments are modified in-place.

        @param matrix: 3x3 array
        @type matrix: array
        @param x_vector: right hand side for X coefficients
        @type x_vector: array
        @param y_vector: right hand side for Y coefficients
        @type y_vector: array
        """

        # positions in the public release (OPO) image, for the matrix
        sum_xo = 0.
        sum_yo = 0.
        sum_xo2 = 0.
        sum_xoyo = 0.
        sum_yo2 = 0.
        # reference and public release images, for the vectors
        sum_xrxo = 0.
        sum_xryo = 0.
        sum_xr = 0.
        sum_yrxo = 0.
        sum_yryo = 0.
        sum_yr = 0.

        for i in range (1, self.total_points + 1):
            if not self.refXY.has_key (i):
                continue                        # a deleted point
            (xr, yr) = self.refXY[i]
            (xo, yo) = self.opoXY[i]
            sum_xo += xo
            sum_yo += yo
            sum_xo2 += xo**2
            sum_xoyo += xo * yo
            sum_yo2 += yo**2
            sum_xrxo += xr * xo
            sum_xryo += xr * yo
            sum_xr += xr
            sum_yrxo += yr * xo
            sum_yryo += yr * yo
            sum_yr += yr

        matrix[0,0] = sum_xo2
        matrix[1,0] = sum_xoyo
        matrix[2,0] = sum_xo
        matrix[0,1] = sum_xoyo
        matrix[1,1] = sum_yo2
        matrix[2,1] = sum_yo
        matrix[0,2] = sum_xo
        matrix[1,2] = sum_yo
        matrix[2,2] = self.npoints

        x_vector[0] = sum_xrxo
        x_vector[1] = sum_xryo
        x_vector[2] = sum_xr

        y_vector[0] = sum_yrxo
        y_vector[1] = sum_yryo
        y_vector[2] = sum_yr

    def computeResiduals (self):
        """Compute residuals of fit in reference image."""

        if self.x_coeff is None:
            return

        sumn = 0
        sumx2 = 0.
        sumy2 = 0.
        for i in range (1, self.total_points + 1):
            if self.refXY.has_key (i):
                (xr, yr) = self.refXY[i]
                (xo, yo) = self.opoXY[i]
                (xr_fit, yr_fit) = self.opoToRef (xo, yo)
                self.residuals[i] = (xr - xr_fit, yr - yr_fit)
                sumn += 1
                sumx2 += (xr - xr_fit)**2
                sumy2 += (yr - yr_fit)**2
        if sumn > 0:
            rms_x = math.sqrt (sumx2 / float (sumn))
            rms_y = math.sqrt (sumy2 / float (sumn))
        else:
            rms_x = 0.
            rms_y = 0.
        self.residuals["rms_x"] = rms_x
        self.residuals["rms_y"] = rms_y

    def opoToRef (self, x, y):
        """Pixels in public release image to pixels in reference image.

        @param x: unbinned pixel coordinates in public release image
            (x is the more rapidly varying axis in the image); the length
            of this array is the number of marked pairs
        @type x: array
        @param y: unbinned pixel coordinates in public release image
            (y is the less rapidly varying axis in the image)
        @type y: array

        @return: a tuple (x_ref, y_ref) of unbinned pixel coordinates in
            the reference image
        @rtype: tuple
        """

        a = self.x_coeff
        b = self.y_coeff

        x_ref = a[0] * x + a[1] * y + a[2]
        y_ref = b[0] * x + b[1] * y + b[2]

        return (x_ref, y_ref)

    def refToOpo (self, x, y):
        """Pixels in reference image to pixels in public release image.

        @param x: unbinned pixel coordinates in reference image
            (x is the more rapidly varying axis in the image); the length
            of this array is the number of marked pairs
        @type x: array
        @param y: unbinned pixel coordinates in reference image
            (y is the less rapidly varying axis in the image)
        @type y: array

        @return: a tuple (x_opo, y_opo) of unbinned pixel coordinates in
            the public-release image
        @rtype: tuple
        """

        if self.x_inverse is None:
            if self.x_coeff is None:
                return (None, None)
            a = self.x_coeff
            b = self.y_coeff
            matrix = N.array (((a[0], a[1]), (b[0], b[1])), dtype=N.float64)
            try:
                inverse = N.linalg.inv (matrix)
            except LinearAlgebraError:
                tkMessageBox.showerror ("singlular matrix",
                    "the fit is singular\nyou must improve the fit first")
                return (None, None)
            self.x_inverse = N.zeros (3, dtype=N.float64)
            self.y_inverse = N.zeros (3, dtype=N.float64)
            self.x_inverse[0] = inverse[0,0]
            self.y_inverse[0] = inverse[1,0]
            self.x_inverse[1] = inverse[0,1]
            self.y_inverse[1] = inverse[1,1]
            self.x_inverse[2] = self.x_coeff[2]
            self.y_inverse[2] = self.y_coeff[2]

        a = self.x_inverse
        b = self.y_inverse
        x_x0 = x - a[2]
        y_y0 = y - b[2]
        x_opo = a[0] * x_x0 + a[1] * y_y0
        y_opo = b[0] * x_x0 + b[1] * y_y0

        return (x_opo, y_opo)

    def xi_eta (self, x, y):
        """Convert from pixel coordinates to tangent-plane coordinates."""

        (xi, eta) = self.ref.wcs.pixelIntermediate ((x, y))

        return (xi, eta)

    def writeWCS (self):
        if self.opo_scale is None:
            tkMessageBox.showerror ("no WCS info",
                "you must determine the WCS info before you can save it")
            return
        default = self.opo.name + ".wcs"
        self.output = tkFileDialog.asksaveasfilename (title="output file",
                parent=self.frame, initialfile=default)
        if not self.output:
            return
        fd = open (self.output, "w")
        self.writeInfo (fd)
        fd.close()
        self.print_wcs = False

    def writeInfo (self, fd):
        fd.write ("public-release image = %s\n" % self.opo.name)
        fd.write ("reference image = %s\n" % self.ref.name)
        fd.write ('WCS_CoordProjection = "TAN"\n')
        # one indexed
        fd.write ("WCS_CoordRefPixel = %.1f, %.1f\n" % \
                  (self.opo.crpix[0] + 1., self.opo.crpix[1] + 1.))
        if self.npoints > 3:
            fd.write ("WCS_CoordRefValue = %11.7f, %10.7f\n" % \
                      (self.opo.crval[0], self.opo.crval[1]))
            fd.write ("WCS_CDMatrix = %.9g, %.9g, %.9g, %.9g\n" % \
                      (self.opo.cd[0,0], self.opo.cd[0,1],
                       self.opo.cd[1,0], self.opo.cd[1,1]))
            fd.write ("scale = %.6g\n" % self.opo_scale)
            fd.write ("orientation = %.1f\n" % self.opo_orientation)
        if self.x_coeff is not None:
            fd.write ("# plate solution:\n")
            fd.write ("x_reference = %.8g * x_pr + %.8g * y_pr + %.8g\n" % \
                      (self.x_coeff[0], self.x_coeff[1], self.x_coeff[2]))
            fd.write ("y_reference = %.8g * x_pr + %.8g * y_pr + %.8g\n" % \
                      (self.y_coeff[0], self.y_coeff[1], self.y_coeff[2]))
        if self.npoints < 1:
             return
        if self.residuals.has_key ("rms_x") and self.residuals["rms_x"] > 0.:
            fd.write ("RMS of fit = %.5g\n" % \
                    math.sqrt (self.residuals["rms_x"]**2 + \
                               self.residuals["rms_y"]**2))
        fd.write ("#  PR image    ref image   evaluate fit  residuals\n")
        #          12345 12345  12345 12345  123456 123456  (123 123)
        for i in range (1, self.total_points + 1):
            if self.refXY.has_key (i):
                (xr, yr) = self.refXY[i]
                (xo, yo) = self.opoXY[i]
                (xr_fit, yr_fit) = self.opoToRef (xo, yo)
                (x_resid, y_resid) = self.residuals[i]
                fd.write (
                "%5.0f %5.0f  %5.0f %5.0f  %6.1f %6.1f  (%3.1f %3.1f)\n" % \
                    (xo, yo, xr, yr, xr_fit, yr_fit, x_resid, y_resid))

    def displayRefImage (self):
        if self.ref.current_data is None:
            return
        pylab.figure (REF_FIGURE)
        xfmt = ticker.FuncFormatter (self.ref.xscale)
        yfmt = ticker.FuncFormatter (self.ref.yscale)
        im_ref = pylab.imshow (self.ref.current_data,
                      aspect="equal", interpolation="nearest", origin="lower",
                      cmap=cm.gray, hold=False,
                      vmin=self.ref.Vmin, vmax=self.ref.Vmax)
        im_ref.axes.xaxis.set_major_formatter (xfmt)
        im_ref.axes.yaxis.set_major_formatter (yfmt)
        im_ref.axes.fmt_xdata = self.ref.fmt_xdata
        im_ref.axes.fmt_ydata = self.ref.fmt_ydata
        pylab.title (os.path.basename (self.ref.name))

    def displayOPOImage (self):
        if self.opo.data is None:
            return
        pylab.figure (OPO_FIGURE)
        xfmt = ticker.FuncFormatter (self.opo.xscale)
        yfmt = ticker.FuncFormatter (self.opo.yscale)
        im_opo = pylab.imshow (self.opo.data,
                      aspect="equal", interpolation="nearest", origin="lower",
                      cmap=cm.gray, hold=False,
                      vmin=self.opo.Vmin, vmax=self.opo.Vmax)
        im_opo.axes.xaxis.set_major_formatter (xfmt)
        im_opo.axes.yaxis.set_major_formatter (yfmt)
        im_opo.axes.fmt_xdata = self.opo.fmt_xdata
        im_opo.axes.fmt_ydata = self.opo.fmt_ydata
        pylab.title (os.path.basename (self.opo.name))

    def read_mouse_xy (self, fig):
        pylab.figure (fig)
        self.done = False
        self.cid = pylab.connect ("button_press_event", self.mouse_xy)
        while not self.done:
            self.frame.update()
            sleep (0.03)

    def mouse_xy (self, event):
        if event.inaxes:
            self.xdata = event.xdata
            self.ydata = event.ydata
            pylab.disconnect (self.cid)
            self.done = True
        if event.button != 1:		# xxx temporary!
            self.addPointsDone = True

class GetRefMinMax (tkSimpleDialog.Dialog):

    def body (self, master):

        global vmin, vmax
        self.vmin = None
        self.vmax = None

        Label (master, text="lower limit:").grid (row=0)
        Label (master, text="upper limit:").grid (row=1)

        self.e1 = Entry (master)
        self.e2 = Entry (master)

        self.e1.grid (row=0, column=1)
        self.e2.grid (row=1, column=1)

        e1 = StringVar()
        e1.set (str (vmin))
        self.e1["textvariable"] = e1
        e2 = StringVar()
        e2.set (str (vmax))
        self.e2["textvariable"] = e2

        return self.e1

    def apply (self):

        self.vmin = string.atof (self.e1.get())
        self.vmax = string.atof (self.e2.get())

class GetOpoMinMax (tkSimpleDialog.Dialog):

    def body (self, master):

        global vomin, vomax
        self.vomin = None
        self.vomax = None

        Label (master, text="lower limit:").grid (row=0)
        Label (master, text="upper limit:").grid (row=1)

        self.e1 = Entry (master)
        self.e2 = Entry (master)

        self.e1.grid (row=0, column=1)
        self.e2.grid (row=1, column=1)

        e1 = StringVar()
        e1.set (str (vmin))
        self.e1["textvariable"] = e1
        e2 = StringVar()
        e2.set (str (vmax))
        self.e2["textvariable"] = e2

        return self.e1

    def apply (self):

        self.vomin = string.atof (self.e1.get())
        self.vomax = string.atof (self.e2.get())



class Image (object):

    def __init__ (self, figure):

        self.figure = figure

        self.name = None                # name of image
        self.data = None                # image data
        self.header = None              # image header

        # length in pixels of each leg of the + sign for marking points
        # (in binned pixels)
        self.marklen = 10       # to be replaced, depending on image size
        self.length_fraction = 0.015    # marklen is this fraction of image

        # half width for finding brightest pixel (in binned pixels)
        self.search = 10        # to be replaced, depending on image size
        self.weight = None      # to be replaced

        # plot limits
        (self.xlow, self.xhigh) = (0., 1.)
        (self.ylow, self.yhigh) = (0., 1.)
        self.Vmin = None                # min data value to display
        self.Vmax = None                # max data value to display

        # WCS info
        # For the reference image, these will be read from the header.
        # For the public-release image, these will be computed.
        # Except for bnaxis, all these are for coordinates in the original,
        # unbinned image.
        # wcs is the celwcs object, while crval, crpix, and cd include only
        # the right ascension and declination information (in that order).
        self.wcs = None
        self.bnaxis = N.zeros (2, dtype=N.int32)   # binned size
        self.naxis = N.zeros (2, dtype=N.int32)
        self.crval = N.zeros (2, dtype=N.float64)
        self.crpix = N.zeros (2, dtype=N.float64)  # zero indexed
        self.cd = N.identity (2, dtype=N.float64)

        # values for the current (possibly rotated) image
        self.rotangle = 0.              # orientation of displayed image
        self.current_data = None                # possibly rotated
        self.current_crpix = self.crpix.copy()
        self.current_cd = self.cd.copy()

        # binned pixel coordinates (i.e. in the internal data arrays)
        # as a function of original pixel coordinates are given by:
        #   x_binned = x_a1 * x + x_a0
        #   y_binned = y_a1 * y + x_a0
        self.block_size = 1
        self.x_a0 = 0.
        self.x_a1 = 1.
        self.y_a0 = 0.
        self.y_a1 = 1.
        # indicates where subsampling was done for block averaging;
        # this will be overridden in refImage and opoImage
        self.bin_locn = LOCN_CORNER

    def setBinning (self, block_size):
        """Compute linear mapping from original to binned pixels.

        The zero point coefficients are for the case that the binning is done
        by boxcar smoothing followed by subsampling by [0,-1,block_size] in
        each axis.
        """

        self.block_size = block_size

        # These coefficients are for a linear mapping from original
        # pixel coordinates to binned pixel coordinates.
        block_size = float (block_size)
        self.x_a1 = 1. / block_size
        self.y_a1 = 1. / block_size

        if self.bin_locn == LOCN_CORNER:
            # appropriate for the case that the binning was done by
            # sampling at [0,0] (lower left corner) in each block
            self.x_a0 = 0.
            self.y_a0 = 0.
        elif self.bin_locn == LOCN_CENTER:
            # appropriate for the case that the binning was done either
            # by sampling at the center or by block averaging
            self.x_a0 = 0.5 * (1./block_size - 1.)
            self.y_a0 = 0.5 * (1./block_size - 1.)

    def toBinned (self, coords):
        """Convert from original pixel coordinates to binned pixels."""

        x = self.x_a1 * coords[0] + self.x_a0
        y = self.y_a1 * coords[1] + self.y_a0
        return (x, y)

    def fromBinned (self, coords_binned):
        """Convert from binned pixel coordinates to original pixels."""

        x = (coords_binned[0] - self.x_a0) / self.x_a1
        y = (coords_binned[1] - self.y_a0) / self.y_a1
        return (x, y)

    # xscale and yscale are used by the "ticker" to display the major axis
    # ticks.  These functions convert from binned pixels to unbinned.
    def xscale (self, x, dummy):
        """Format the X pixel coordinate for the axis ticks."""
        return "%4d" % int ((x - self.x_a0) / self.x_a1)

    def yscale (self, y, dummy):
        """Format the Y pixel coordinate for the axis ticks."""
        return "%4d" % int ((y - self.y_a0) / self.y_a1)

    # fmt_xdata and fmt_ydata are used for displaying pixel coordinates
    # on the menu bar.  These convert from binned pixels to unbinned, and
    # the "%d" format shows the pixel numbers to full precision.
    def fmt_xdata (self, x):
        """Convert to unbinned pixels, and format for display."""
        return "%d" % int ((x - self.x_a0) / self.x_a1)

    def fmt_ydata (self, y):
        """Convert to unbinned pixels, and format for display."""
        return "%d" % int ((y - self.y_a0) / self.y_a1)

    def setSearchParameters (self):
        """Assign values for the search width and weight array."""

        # the half width of the search region
        self.search = int (round (self.marklen / 2.))

        # quadratic weighting
        shape = list (self.data.shape)          # may be rank 2 or rank 3
        shape[0] = 2 * self.search + 1  # this is the size of the search region
        shape[1] = 2 * self.search + 1
        nelem = 1
        for n in shape:
            nelem *= n
        self.weight = N.ones (shape=shape, dtype=N.float32)
        srch = float (self.search)
        for j in range (2*self.search+1):
            for i in range (2*self.search+1):
                wgt = (srch + 1.)**2 - ((i - srch)**2 + (j - srch)**2)
                self.weight[j,i] = max (wgt, 0.)

    def drawPlus (self, x0=0., y0=0.,
                  index=0, color="w", lenfactor=1,
                  draw_it=True, rescale=False):
        """Draw a plus sign with a central gap at (x0, y0).

        The input pixel coordinates will be converted to binned pixel
        coordinates, and the plus sign will actually be drawn at the
        latter location.

        @param x0: unbinned pixel coordinate (more rapidly varying axis)
        @type x0: float
        @param y0: unbinned pixel coordinate (less rapidly varying axis)
        @type y0: float
        @param index: if this integer is greater than zero, also draw this
            number
        @type index: int
        @param color: the color to make the plus sign and optional number
        @type color: string
        @param lenfactor: multiply self.marklen by this to get length to plot
        @type lenfactor: int
        @param draw_it: if true, redraw the plot
        @type draw_it: boolean
        @param rescale: if true, reset xlim and ylim for the plot
        @type rescale: boolean
        """

        pylab.figure (self.figure)

        length = int (lenfactor * self.marklen)
        gap = 3 * length // 10

        # convert to binned pixels
        (x0_b, y0_b) = self.toBinned ((x0, y0))
        x0_b += 0.5
        y0_b += 0.5

        x = N.zeros (8, dtype=N.float64)
        y = N.zeros (8, dtype=N.float64)
        # horizontal part
        y[0:4] = y0_b
        x[0] = x0_b - length
        x[1] = x0_b - gap
        x[2] = x0_b + gap
        x[3] = x0_b + length
        # vertical part
        x[4:8] = x0_b
        y[4] = y0_b - length
        y[5] = y0_b - gap
        y[6] = y0_b + gap
        y[7] = y0_b + length
        ax = pylab.gca()
        ax.plot (x[0:2], y[0:2], color=color)
        ax.plot (x[2:4], y[2:4], color=color)
        ax.plot (x[4:6], y[4:6], color=color)
        ax.plot (x[6:8], y[6:8], color=color)

        if index > 0:
            ax.text (x0_b + length + gap, y0_b, str (index), color=color,
                     horizontalalignment="center",
                     verticalalignment="center")
        if draw_it:
            pylab.draw()

        if rescale:
            pylab.xlim (self.xlow, self.xhigh)
            pylab.ylim (self.ylow, self.yhigh)

def computewcs():

    root = Tk()
    wcs = ComputeWCS (root)
    root.mainloop()

if __name__ == "__main__":

    computewcs()
