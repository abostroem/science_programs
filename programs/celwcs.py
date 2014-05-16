#! /usr/bin/env python

import types
import math
import re
import copy
import numpy as N
import pyfits

__version__ = "2007 September 5"

DEGREEStoRADIANS = math.pi / 180.       # convert degrees to radians
TINY = 1.e-10                           # nearly zero

# WCS keywords have different names for these three different ways
# of storing image data.  The name gives the term by which these
# are identified in the FITS documents.
PRIMARY_ARRAY  = 1      # a primary HDU or an IMAGE HDU
BINTABLE_ARRAY = 2      # an array in one row of a column of arrays
PIXEL_LIST     = 3      # a tabulated list of pixels, using two columns

# Distinguish between zenithal projections and others, because the
# default values for the longitude and latitude of the fiducial point
# are not the same.
ZENITHAL_PROJ = ["AZP", "SZP", "TAN", "STG", "SIN", "ARC", "ZPN", "ZEA", "AIR"]
CYLINDRICAL_PROJ = ["CYP", "CEA", "CAR", "MER", "SFL", "PAR", "MOL", "AIT"]
CONIC_PROJ = ["COP", "COE", "COD", "COO"]
OTHER_PROJ = ["BON", "PCO", "TSC", "CSC", "QSC", "HPX"]
SPHERICAL_MAP_PROJ = ZENITHAL_PROJ + CYLINDRICAL_PROJ + CONIC_PROJ + OTHER_PROJ

def readCoords (input, col=(0,1), extension=1):
    """Read coordinates from a FITS table.

    For a table containing these three coordinate vectors:
        x1  y1
        x2  y2
        x3  y3
    the output would be:
       [[x1  x2  x3]
        [y1  y2  y3]]
    If the table contained only the first row, the output would be:
       [x1  y1]

    @param input: name of a FITS file containing a bintable (or table).
    @type input: string
    @param col: the names or numbers (zero indexed) of columns in the table
        containing the coordinates; the default is to read the first two
        columns as X and Y.
    @type col: tuple
    @param extension: the number or EXTNAME of the extension containing
        the table.
    @type extension: int or string

    @return: an array of coordinates read from the table, in a format suitable
        for use as input to either frompixel() or topixel(); if there is only
        one row in the table, the function value will be a single vector
    @rtype: ndarray
    """

    if type (extension) == types.IntType and extension < 1:
        raise ValueError, "'extension' must be greater than zero"

    fd = pyfits.open (input)
    if type (extension) == types.IntType and len (fd) < extension+1:
        fd.close()
        raise ValueError, "extension %d does not exist in input file" % \
                          extension
    if fd[extension].data is None or len (fd[extension].data) < 1:
        fd.close()
        errmess = "no data in extension " + str (extension)
        raise ValueError, errmess

    nrows = len (fd[extension].data)

    coldata = []
    for column in col:
        if nrows == 1:
            coldata.append (fd[extension].data.field (column)[0])
        else:
            coldata.append (fd[extension].data.field (column))

    return N.array (coldata)

##############################################################################
class WCSbase(object):
    """Base object for WCS objects.

    Contains basic functionality such as header parsing."""
    def __init__ (self, hdr, alternate=None, column=None):

        # this is a list of all WCS keywords in hdr
        self.keywords = []

        # set keyword_translate attribute
        self._setKeywordDict()

        # set attributes for names of spherical projection functions
        # and (eventually) distortion functions
        self._setFunctionAttributes()

        # copy arguments to attributes, with some checking and interpretation
        self._saveArguments (hdr, alternate, column)

        # find the dimension of the WCS from a keyword
        self._findWcsAxes (hdr)

        # create attribute arrays with the right size and default values
        self._initializeWcsAttributes()

        # copy keyword values to attributes
        self._extractWcsKeywords (hdr)

    def _setKeywordDict (self):
        """Set the mapping from header keyword to attribute name."""

        self.keyword_translate = {
            "WCSAXES": "wcsaxes",       # primary array keywords
            "CTYPE":   "ctype",
            "CUNIT":   "cunit",
            "CRVAL":   "crval",
            "CDELT":   "cdelt",
            "CRPIX":   "crpix",
            "PC":      "pc",
            "CD":      "cd",
            "PV":      "pv",            # also bintable array
            "PS":      "ps",            # also bintable array
            "WCSNAME": "wcsname",
            "CRDER":   "crder",
            "CSYER":   "csyer",
            "CROTA":   "crota",
            "LONPOLE": "lonpole",
            "LATPOLE": "latpole",
            "EQUINOX": "equinox",
            "MJD-OBS": "mjd_obs",
            "RADESYS": "radesys",
            "WCAX":    "wcsaxes",       # bintable array
            "CTYP":    "ctype",
            "CUNI":    "cunit",
            "CRVL":    "crval",
            "CDLT":    "cdelt",
            "CRPX":    "crpix",
            "CTY":     "ctype",
            "CUN":     "cunit",
            "CRV":     "crval",
            "CDE":     "cdelt",
            "CRP":     "crpix",
            "V":       "pv",
            "S":       "ps",
            "CRD":     "crder",
            "CSY":     "csyer",
            "CROT":    "crota",
            "WCST":    "wcst",
            "WCSX":    "wcsx",
            "LONP":    "lonpole",       # both bintable array and pixel list
            "LATP":    "latpole",
            "EQUI":    "equinox",
            "MJDOB":   "mjd_obs",
            "RADE":    "radesys",
            "WCSN":    "wcsname",
            "TCTYP":   "ctype",         # pixel list
            "TCUNI":   "cunit",
            "TCRVL":   "crval",
            "TCDLT":   "cdelt",
            "TCRPX":   "crpix",
            "TCTY":    "ctype",
            "TCUN":    "cunit",
            "TCRV":    "crval",
            "TCDE":    "cdelt",
            "TCRP":    "crpix",
            "TP":      "pc",
            "TPC":     "pc",
            "TC":      "cd",
            "TCD":     "cd",
            "TV":      "pv",
            "TPV":     "pv",
            "TS":      "ps",
            "TPS":     "ps",
            "TCRD":    "crder",
            "TCSY":    "csyer",
            "TCROT":   "crota"}

    def _setFunctionAttributes (self):
        """Assign the attributes for spherical projection functions."""

        # These will be set to one pair of the functions listed below.
        self._projection_fcn = None
        self._inverse_fcn = None

        self.TAN_projection = self._tanProj
        self.TAN_inverse = self._tanInvProj

        # xxx not implemented yet
        self.AZP_projection = self._azpProj
        self.AZP_inverse = self._azpInvProj

        self.STG_projection = self._stgProj
        self.STG_inverse = self._stgInvProj

        self.SIN_projection = self._sinProj
        self.SIN_inverse = self._sinInvProj

        self.ARC_projection = self._arcProj
        self.ARC_inverse = self._arcInvProj

        self.SFL_projection = self._sflProj
        self.SFL_inverse = self._sflInvProj

        self.PAR_projection = self._parProj
        self.PAR_inverse = self._parInvProj

        self.MOL_projection = self._molProj
        self.MOL_inverse = self._molInvProj

        self.AIT_projection = self._aitProj
        self.AIT_inverse = self._aitInvProj

        self.MER_projection = self._merProj
        self.MER_inverse = self._merInvProj

        self.COP_projection = self._copProj
        self.COP_inverse = self._copInvProj

        self.COE_projection = self._coeProj
        self.COE_inverse = self._coeInvProj

        self.BON_projection = self._bonProj
        self.BON_inverse = self._bonInvProj

    def _tanProj (self, phi, theta):
        """gnomonic (tangent) projection"""

        r = N.tan (self.lat_0 - theta)
        x = r * N.sin (phi)
        y = -r * N.cos (phi)

        return (x, y)

    def _tanInvProj (self, x, y):
        """inverse gnomonic (tangent) projection"""

        r = N.sqrt (x**2 + y**2)
        theta = self.lat_0 - N.arctan (r)
        phi = N.arctan2 (x, -y)

        return (phi, theta)

    def _azpProj (self, phi, theta):
        """zenithal perspective projection"""

        # xxx not implemented yet
        print "\nERROR: AZP not yet implemented"
        r = N.tan (self.lat_0 - theta)
        x = r * N.sin (phi)
        y = -r * N.cos (phi)

        return (x, y)

    def _azpInvProj (self, x, y):
        """inverse zenithal perspective projection"""

        # xxx not implemented yet
        print "\nERROR: Inverse AZP not yet implemented"
        r = N.sqrt (x**2 + y**2)
        theta = self.lat_0 - N.arctan (r)
        phi = N.arctan2 (x, -y)

        return (phi, theta)

    def _stgProj (self, phi, theta):
        """stereographic projection"""

        r = N.tan ((self.lat_0 - theta) / 2.)
        x = r * N.sin (phi)
        y = -r * N.cos (phi)

        return (x, y)

    def _stgInvProj (self, x, y):
        """inverse stereographic projection"""

        r = N.sqrt (x**2 + y**2)
        theta = self.lat_0 - 2. * N.arctan (r)
        phi = N.arctan2 (x, -y)

        return (phi, theta)

    def _sinProj (self, phi, theta):
        """sine projection"""

        r = N.cos (theta)
        x = r * N.sin (phi)
        y = -r * N.cos (phi)

        return (x, y)

    def _sinInvProj (self, x, y):
        """inverse sine projection"""

        r = N.sqrt (x**2 + y**2)
        theta = N.arccos (r)
        phi = N.arctan2 (x, -y)

        return (phi, theta)

    def _arcProj (self, phi, theta):
        """zenithal equidistant projection"""

        r = self.lat_0 - theta
        x = r * N.sin (phi)
        y = -r * N.cos (phi)

        return (x, y)

    def _arcInvProj (self, x, y):
        """inverse zenithal equidistant projection"""

        r = N.sqrt (x**2 + y**2)
        theta = self.lat_0 - r
        phi = N.arctan2 (x, -y)

        return (phi, theta)

    def _sflProj (self, phi, theta):
        """Sanson-Flamsteed projection"""

        x = phi * N.cos (theta)
        y = theta

        return (x, y)

    def _sflInvProj (self, x, y):
        """inverse Sanson-Flamsteed projection"""

        theta = y
        phi = x / N.cos (y)

        return (phi, theta)

    def _parProj (self, phi, theta):
        """parabolic projection"""

        x = phi * (2. * N.cos (2. * theta / 3.) - 1.)
        y = N.pi * N.sin (theta / 3.)

        return (x, y)

    def _parInvProj (self, x, y):
        """inverse parabolic projection"""

        theta = 3. * N.arcsin (y / N.pi)
        phi = x / (1. - 4. * (y / N.pi)**2)

        return (phi, theta)

    def _molProj (self, phi, theta):
        """Mollweide's projection"""

        gamma = self.molGamma (theta, nloops=10)
        x = phi * (2. * N.sqrt (2.) / N.pi) * N.cos (gamma)
        y = N.sqrt (2.) * N.sin (gamma)

        return (x, y)

    def _molInvProj (self, x, y):
        """inverse Mollweide's projection"""

        theta = N.arcsin (2. * N.arcsin (y / N.sqrt (2.)) / N.pi + \
                           y * N.sqrt (2. - y**2) / N.pi)
        phi = N.pi * x / (2. * N.sqrt (2. - y**2))

        return (phi, theta)

    def _molGamma (self, theta, nloops=10):
        """Compute gamma for Mollweide's projection

        Given theta, we need to solve the following for gamma:
            sin (theta) = (2*gamma + sin (2*gamma)) / pi
        """

        gamma = theta * N.pi / 4.       # first approximation

        for i in range (nloops):
            test_theta = N.arcsin ((2. * gamma + N.sin (2. * gamma)) / N.pi)
            slope = (N.pi / 2.) * N.cos (test_theta) / \
                    (1. + N.cos (2. * gamma))
            gamma += (theta - test_theta) * slope

        return gamma

    def _aitProj (self, phi, theta):
        """Hammer-Aitoff projection"""

        gamma = N.sqrt (2. / (1. + N.cos (theta) * N.cos (phi / 2.)))
        x = 2. * gamma * N.cos (theta) * N.sin (phi / 2.)
        y = gamma * N.sin (theta)

        return (x, y)

    def _aitInvProj (self, x, y):
        """inverse Hammer-Aitoff projection"""

        z = N.sqrt (1. - (x/4.)**2 - (y/2.)**2)
        theta = N.arcsin (y * z)
        phi = 2. * N.arctan2 (x * z/2., 2.*z**2 - 1.)

        return (phi, theta)

    def _merProj (self, phi, theta):
        """Mercator projection"""

        x = phi
        y = N.log (N.tan ((N.pi/2. + theta) / 2.))

        return (x, y)

    def _merInvProj (self, x, y):
        """inverse Mercator projection"""

        theta = 2. * N.arctan (N.exp (y)) - N.pi / 2.
        phi = x

        return (phi, theta)

    def _copProj (self, phi, theta):
        """conic perspective projection"""

        (C, cot_theta_a, cos_eta, y0) = self.copInit()
        r = cos_eta * (cot_theta_a - N.tan (theta - self.theta_a))

        x = r * N.sin (C * phi)
        y = -r * N.cos (C * phi) + y0

        return (x, y)

    def _copInvProj (self, x, y):
        """inverse conic perspective projection"""

        (C, cot_theta_a, cos_eta, y0) = self.copInit()
        r = N.sign (self.theta_a) * N.sqrt (x**2 + (y0 - y)**2)

        theta = self.theta_a + N.arctan (cot_theta_a - r / cos_eta)
        phi = N.arctan2 (x / r, (y0 - y) / r) / C

        return (phi, theta)

    def _copInit (self):
        """Compute some parameters for conic perspective projection"""

        if self.theta_a is None:
            raise RuntimeError, "theta_a is undefined " \
                  "(use the appropriate PV keyword)"
        if self.theta_a == 0. or abs (self.theta_a) == math.pi / 2.:
            raise ValueError, "theta_a = %.10g (from PV keyword) " \
                  "is not supported for conic projection" % self.theta_a
        C = math.sin (self.theta_a)
        cot_theta_a = 1. / math.tan (self.theta_a)
        cos_eta = math.cos (self.eta)
        y0 = cos_eta * cot_theta_a

        return (C, cot_theta_a, cos_eta, y0)

    def _coeProj (self, phi, theta):
        """conic equal area projection"""

        (C, sin_theta1, sin_theta2, gamma, y0) = self.coeInit()
        r = (2. / gamma) * N.sqrt (1. + \
                 sin_theta1 * sin_theta2 - gamma * N.sin (theta))

        x = r * N.sin (C * phi)
        y = -r * N.cos (C * phi) + y0

        return (x, y)

    def _coeInvProj (self, x, y):
        """inverse conic equal area projection"""

        (C, sin_theta1, sin_theta2, gamma, y0) = self.coeInit()
        r = N.sign (self.theta_a) * N.sqrt (x**2 + (y0 - y)**2)

        theta = N.arcsin (1. / gamma + sin_theta1 * sin_theta2 / gamma - \
                           gamma * (r/2.)**2)
        phi = N.arctan2 (x / r, (y0 - y) / r) / C

        return (phi, theta)

    def _coeInit (self):
        """Compute some parameters for conic equal area projection"""

        if self.theta_a is None:
            raise RuntimeError, "theta_a is undefined " \
                  "(use the appropriate PV keyword)"

        theta1 = self.theta_a - self.eta
        theta2 = self.theta_a + self.eta
        sin_theta1 = math.sin (theta1)
        sin_theta2 = math.sin (theta2)
        gamma = sin_theta1 + sin_theta2
        C = gamma / 2.
        y0 = (2. / gamma) * math.sqrt (1. + \
                sin_theta1 * sin_theta2 - gamma * math.sin (self.theta_a))

        return (C, sin_theta1, sin_theta2, gamma, y0)

    def _bonProj (self, phi, theta):
        """Bonne's equal area projection"""

        y0 = self.theta_a + 1. / math.tan (self.theta_a)
        r = y0 - theta
        A = (phi / r) * N.cos (theta)

        x = r * N.sin (A)
        y = -r * N.cos (A) + y0

        return (x, y)

    def _bonInvProj (self, x, y):
        """inverse Bonne's equal area projection"""

        y0 = self.theta_a + 1. / math.tan (self.theta_a)
        r = N.sign (self.theta_a) * N.sqrt (x**2 + (y0 - y)**2)
        A = N.arctan2 (x / r, (y0 - y) / r)

        theta = y0 - r
        phi = A * r / N.cos (theta)

        return (phi, theta)

    def _saveArguments (self, hdr, alternate, column):
        """Check arguments, and copy to attributes.

        @param hdr: a set of FITS keywords and values (values are case
            sensitive, but keyword names are not)
        @type hdr: pyfits.Header object or dictionary
        @param alternate: a single letter, specifying an alternate WCS
            (default is blank, for primary WCS)
        @type alternate: string (one letter) or None
        @param column: for an image array, this should be None;
            for an array in a table column, this should be a single column
            number (zero indexed);
            for a tabulated list of pixels in two columns, this should be a
            list giving the column numbers for right ascension (longitude)
            and declination (latitude)
        @type column: None, int, or list of two ints
        """

        if not isinstance (hdr, (pyfits.Header, types.DictionaryType)):
            raise TypeError, \
        "'hdr' must be a Header object or a dictionary of keywords and values"

        self.hdr = hdr

        if alternate is None or len (alternate) < 1 or alternate == " ":
            self.alternate = ""
        elif len (alternate) > 1:
            raise ValueError, "'alternate' must be a single character or blank"
        else:
            self.alternate = alternate.upper()

        if column is None:
            self.image_rep = PRIMARY_ARRAY
            self.column = None
            self.tab_columns = None
        else:
            if isinstance (column, (types.ListType, types.TupleType)):
                if len (column) < 2:
                    raise ValueError, \
        " when 'column' is a list it must contain two or more column numbers"
                elif column[0] < 0 or column[1] < 0:
                    raise ValueError, \
                        "the column numbers in 'column' must be non-negative"
                self.image_rep = PIXEL_LIST
                self.column = None
                self.tab_columns = list (column)
            else:
                if column < 0:
                    raise ValueError, \
                        "column number 'column' must be non-negative"
                self.image_rep = BINTABLE_ARRAY
                self.column = column
                self.tab_columns = None

    def _findWcsAxes (self, hdr):
        """Find the dimension of the WCS system.

        @param hdr: a set of FITS keywords and values (values are case
            sensitive, but keyword names are not)
        @type hdr: pyfits.Header object or dictionary

        The dimension will be taken from keyword WCSAXES if it is present, or
        NAXIS, or the dimension will be set to 2 if neither naxis nor wcsaxes
        is found in hdr.
        For a pixel list, the user is expected to supply a list of column
        numbers, one for each axis, so the length of that list is the
        dimension of the WCS.
        """

        self.naxis = 2                  # default

        items = hdr.items()

        if self.image_rep == PIXEL_LIST:
            self.wcsaxes = len (self.tab_columns)
            for (keyword, value) in items:
                if keyword.upper() == "NAXIS":
                    self.naxis = value
                    break
            if self.wcsaxes < self.naxis:
                self.wcsaxes = self.naxis
            return

        if self.image_rep == PRIMARY_ARRAY:
            r = re.compile (r"""(?x)
                (?P<name>WCSAXES)
                (?P<alt>[A-Z]?$)
            """)
        elif self.image_rep == BINTABLE_ARRAY:
            r = re.compile (r"""(?x)
                (?P<name>WCAX)
                (?P<col>\d+)
                (?P<alt>[A-Z]?$)
            """)

        self.wcsaxes = None
        for (keyword, value) in items:
            keyword = keyword.upper()
            if keyword == "NAXIS":
                self.naxis = value
                continue
            m = r.match (keyword)
            if m and m.group ("alt") == self.alternate:
                if self.image_rep == BINTABLE_ARRAY:
                    if m.group ("col") == "" or \
                            int (m.group ("col")) != (self.column + 1):
                        continue
                self.wcsaxes = value
                self.keywords.append (keyword)
                break

        if self.wcsaxes is None:
            self.wcsaxes = self.naxis

        if self.wcsaxes < self.naxis:
            print "Warning:  wcsaxes = %d is smaller than naxis = %d;" % \
                  (self.wcsaxes, self.naxis)
            print "wcsaxes will be set to naxis."
            self.wcsaxes = self.naxis

    def _extractWcsKeywords (self, hdr):
        """Copy WCS keyword values to attributes.

        @param hdr: a set of FITS keywords and values (values are case
            sensitive, but keyword names are not)
        @type hdr: pyfits.Header object or dictionary
        """

        if self.image_rep == PRIMARY_ARRAY:
            self._imageArrayKeywords (hdr)
        elif self.image_rep == BINTABLE_ARRAY:
            self._bintableArrayKeywords (hdr)
        elif self.image_rep == PIXEL_LIST:
            self._pixelListKeywords (hdr)

    def _checkForThese (self, key, keywords):
        """Set an attribute (a flag) if key is in keywords."""

        if key in keywords:
            attrib = "got_" + key
            self.__setattr__ (attrib, True)

    def _imageArrayKeywords (self, hdr):
        """Copy WCS keyword values to attributes.

        @param hdr: a set of FITS keywords and values (values are case
            sensitive, but keyword names are not)
        @type hdr: pyfits.Header object or dictionary

        This version is for primary array keywords.
        """

        r1 = re.compile (r"""(?x)
             (?P<name>WCSNAME|LONPOLE|LATPOLE|EQUINOX|RADESYS)
             (?P<alt>[A-Z]?$)                   # optional alternate WCS
             """)

        r2 = re.compile (r"""(?x)
             (?P<name>CTYPE|CUNIT|CRVAL|CDELT|CRPIX|
                      CRDER|CSYER|CROTA)
             (?P<index>\d+)                     # axis number
             (?P<alt>[A-Z]?$)
             """)

        r3 = re.compile (r"""(?x)
             (?P<name>PC|CD)
             (?P<index1>\d+)_(?P<index2>\d+)    # axis numbers
             (?P<alt>[A-Z]?$)
             """)

        r4 = re.compile (r"""(?x)
             (?P<name>PV|PS)
             (?P<index>\d+)_(?P<par>\d+)        # axis and parameter numbers
             (?P<alt>[A-Z]?$)
             """)

        items = hdr.items()

        for (keyword, value) in items:

            keyword = keyword.upper()
            if keyword == "HISTORY" or keyword == "COMMENT" or keyword == "":
                continue

            if keyword == "MJD-OBS":
                self.mjd_obs = value
                self.keywords.append (keyword)
                continue

            m = r1.match (keyword)
            if m:
                key = self.keyword_translate[m.group ("name")]
                if m.group ("alt") == self.alternate:
                    self.__setattr__ (key, value)
                    self.keywords.append (keyword)
                continue

            m = r2.match (keyword)
            if m and m.group ("alt") == self.alternate:
                key = self.keyword_translate[m.group ("name")]
                if key == "crpix":
                    value -= 1.                 # zero indexed
                a = self.__getattribute__ (key)
                # WTB: hack for legacy CROTA keyword w/o index
                if((m.group("index")=='') and (key=="crota")):
                    print "Warning:  keyword CROTA (with no index) was found"
                    a[0],a[1]=0.0,value # assign CROTA as CROTA2
                    self.__setattr__(key,a)
                else:
                    i = int (m.group ("index")) - 1
                    a[i] = value
                    self.__setattr__ (key, a)
                self.keywords.append (keyword)
                self._checkForThese (key, ("cdelt", "crota"))
                continue

            m = r3.match (keyword)
            if m and m.group ("alt") == self.alternate:
                key = self.keyword_translate[m.group ("name")]
                a = self.__getattribute__ (key)
                i = int (m.group ("index1")) - 1
                j = int (m.group ("index2")) - 1
                a[i,j] = value
                self.__setattr__ (key, a)
                self.keywords.append (keyword)
                self._checkForThese (key, ("cd", "pc"))
                continue

            m = r4.match (keyword)
            if m and m.group ("alt") == self.alternate:
                key = self.keyword_translate[m.group ("name")]
                i = int (m.group ("index")) - 1
                par = int (m.group ("par"))     # already zero indexed
                if key == "pv":
                    self.npv += 1
                    self.pv.append ((i, par, value))
                    self.keywords.append (keyword)
                elif key == "ps":
                    self.nps += 1
                    self.ps.append ((i, par, value))
                    self.keywords.append (keyword)
                continue

        # Check for RADECSYS if RADESYS not found.
        if self.radesys is None and hdr.has_key ("radecsys"):
            self.radesys = hdr["radecsys"]
            self.keywords.append (keyword)

    def _bintableArrayKeywords (self, hdr):
        """Copy WCS keyword values to attributes.

        @param hdr: a set of FITS keywords and values (values are case
            sensitive, but keyword names are not)
        @type hdr: pyfits.Header object or dictionary

        This version is for bintable array keywords.
        """

        r0 = re.compile (r"""(?x)
            (?P<name>MJDOB)
            (?P<col>\d+$)                       # column number
        """)

        r1 = re.compile (r"""(?x)
            (?P<name>LONP|LATP|EQUI|RADE|WCSN)
            (?P<col>\d+)                        # column number
            (?P<alt>[A-Z]?$)
        """)

        r2 = re.compile (r"""(?x)
            (?P<index>\d)                       # axis number
            (?P<name>CTYP|CUNI|CRVL|CDLT|CRPX|
                     CTY|CUN|CRV|CDE|CRP|
                     CRD|CSY|WCST|WCSX|CROT)
            (?P<col>\d+)
            (?P<alt>[A-Z]?$)
        """)

        r3 = re.compile (r"""(?x)
            (?P<index1>\d)                      # first axis number
            (?P<index2>\d)                      # second axis number
            (?P<name>PC|CD)
            (?P<col>\d+)
            (?P<alt>[A-Z]?$)
        """)

        r4 = re.compile (r"""(?x)
            (?P<index>\d)                       # axis number
            (?P<name>V|PV|S|PS)
            (?P<col>\d+)_(?P<par>\d+)           # column and parameter numbers
            (?P<alt>[A-Z]?$)
        """)

        items = hdr.items()

        for (keyword, value) in items:

            keyword = keyword.upper()
            if keyword == "HISTORY" or keyword == "COMMENT" or keyword == "":
                continue

            m = r0.match (keyword)      # this is for MJDOBn (MJD-OBS)
            if m:
                if int (m.group ("col")) != (self.column + 1):
                    continue
                self.mjd_obs = value
                self.keywords.append (keyword)
                continue

            m = r1.match (keyword)
            if m:
                if int (m.group ("col")) != (self.column + 1):
                    continue
                key = self.keyword_translate[m.group ("name")]
                if m.group ("alt") == self.alternate:
                    self.__setattr__ (key, value)
                    self.keywords.append (keyword)
                continue

            m = r2.match (keyword)
            if m:
                if int (m.group ("col")) != (self.column + 1) or \
                        m.group ("alt") != self.alternate:
                    continue
                key = self.keyword_translate[m.group ("name")]
                if key == "crpix":
                    value -= 1.                 # zero indexed
                a = self.__getattribute__ (key)
                i = int (m.group ("index")) - 1
                a[i] = value
                self.__setattr__ (key, a)
                self.keywords.append (keyword)
                self._checkForThese (key, ("cdelt", "crota"))
                continue

            m = r3.match (keyword)
            if m:
                if int (m.group ("col")) != (self.column + 1) or \
                        m.group ("alt") != self.alternate:
                    continue
                key = self.keyword_translate[m.group ("name")]
                a = self.__getattribute__ (key)
                i = int (m.group ("index1")) - 1
                j = int (m.group ("index2")) - 1
                a[i,j] = value
                self.__setattr__ (key, a)
                self.keywords.append (keyword)
                self._checkForThese (key, ("cd", "pc"))
                continue

            m = r4.match (keyword)
            if m:
                if int (m.group ("col")) != (self.column + 1) or \
                        m.group ("alt") != self.alternate:
                    continue
                key = self.keyword_translate[m.group ("name")]
                i = int (m.group ("index")) - 1
                par = int (m.group ("par"))     # already zero indexed
                if key == "pv":
                    self.npv += 1
                    self.pv.append ((i, par, value))
                    self.keywords.append (keyword)
                elif key == "ps":
                    self.nps += 1
                    self.ps.append ((i, par, value))
                    self.keywords.append (keyword)
                continue

    def _pixelListKeywords (self, hdr):
        """Copy WCS keyword values to attributes.

        @param hdr: a set of FITS keywords and values (values are case
            sensitive, but keyword names are not)
        @type hdr: pyfits.Header object or dictionary

        This version is for pixel list keywords.
        """

        r0 = re.compile (r"""(?x)
            (?P<name>MJDOB)
            (?P<col>\d+$)                       # column number
        """)

        r1 = re.compile (r"""(?x)
            (?P<name>LONP|LATP|EQUI|RADE|WCSN)
            (?P<col>\d+)                        # column number
            (?P<alt>[A-Z]?$)
        """)

        r2 = re.compile (r"""(?x)
            (?P<name>TCTYP|TCUNI|TCRVL|TCDLT|TCRPX|
                     TCTY|TCUN|TCRV|TCDE|TCRP|
                     TCRD|TCSY|TCROT)
            (?P<col>\d+)                        # column number
            (?P<alt>[A-Z]?$)
        """)

        r3 = re.compile (r"""(?x)
            (?P<name>TP|TPC|TC|TCD)
            (?P<col1>\d+)_(?P<col2>\d+)         # column numbers
            (?P<alt>[A-Z]?$)
        """)

        r4 = re.compile (r"""(?x)
            (?P<name>TV|TPV|TS|TPS)
            (?P<col>\d+)_(?P<par>\d+)           # column and parameter numbers
            (?P<alt>[A-Z]?$)
        """)

        items = hdr.items()

        for (keyword, value) in items:

            keyword = keyword.upper()
            if keyword == "HISTORY" or keyword == "COMMENT" or keyword == "":
                continue

            m = r0.match (keyword)      # this is for MJDOBn (MJD-OBS)
            if m:
                if (int (m.group ("col")) - 1) not in self.tab_columns:
                    continue
                self.mjd_obs = value
                self.keywords.append (keyword)
                continue

            m = r1.match (keyword)
            if m:
                if (int (m.group ("col")) - 1) not in self.tab_columns:
                    continue
                key = self.keyword_translate[m.group ("name")]
                if m.group ("alt") == self.alternate:
                    self.__setattr__ (key, value)
                    self.keywords.append (keyword)
                continue

            m = r2.match (keyword)
            if m and m.group ("alt") == self.alternate:
                col = int (m.group ("col")) - 1
                if self.tab_columns.count (col) < 1:
                    continue
                i = self.tab_columns.index (col)
                key = self.keyword_translate[m.group ("name")]
                if key == "crpix":
                    value -= 1.                 # zero indexed
                a = self.__getattribute__ (key)
                a[i] = value
                self.__setattr__ (key, a)
                self.keywords.append (keyword)
                self._checkForThese (key, ("cdelt", "crota"))
                continue

            m = r3.match (keyword)
            if m and m.group ("alt") == self.alternate:
                col1 = int (m.group ("col1")) - 1
                col2 = int (m.group ("col2")) - 1
                if self.tab_columns.count (col1) < 1 or \
                   self.tab_columns.count (col2) < 1 :
                    continue
                i = self.tab_columns.index (col1)
                j = self.tab_columns.index (col2)
                key = self.keyword_translate[m.group ("name")]
                a = self.__getattribute__ (key)
                a[i,j] = value
                self.__setattr__ (key, a)
                self.keywords.append (keyword)
                self._checkForThese (key, ("cd", "pc"))
                continue

            m = r4.match (keyword)
            if m and m.group ("alt") == self.alternate:
                col = int (m.group ("col")) - 1
                if self.tab_columns.count (col) < 1:
                    continue
                i = self.tab_columns.index (col)
                key = self.keyword_translate[m.group ("name")]
                par = int (m.group ("par"))     # already zero indexed
                if key == "pv":
                    self.npv += 1
                    self.pv.append ((i, par, value))
                    self.keywords.append (keyword)
                elif key == "ps":
                    self.nps += 1
                    self.ps.append ((i, par, value))
                    self.keywords.append (keyword)
                continue

    def _padPixel (self, pixel):
        """If len (pixel) < wcsaxes, pad the pixel coordinates with zeros.

        @param pixel: 
        @type pixel: ndarray

        @return: 
        @rtype: ndarray
        """

        length = len (pixel)
        if length != self.naxis and length != self.wcsaxes:
            raise ValueError, \
                "The pixel coordinate array must be of length NAXIS or WCSAXES"

        if length == self.wcsaxes:
            return pixel

        shape = pixel.shape
        if len (shape) != 2:
            raise ValueError, \
                "The pixel coordinate array must be either 1-D or 2D"
        padshape = list (shape)
        padshape[0] = self.wcsaxes - length
        pad = N.zeros (padshape, dtype=N.float64)

        return N.concatenate ((pixel, pad), axis=0)

    def _lin_p2i (self, pixel):
        """Do the linear projection from pixels.

        @param pixel: one or more pixel position vectors in the form
            [[x]   or   [[x0  x1  x2]
             [y]]        [y0  y1  y2]]
        @type pixel: ndarray (rank 2)

        @return: intermediate pixel coordinates (i.e. excluding the scale
            factor to physical units):
            output = PC_matrix * (pixel - crpix)
        @rtype: ndarray (rank 2)
        """

        intermediate = pixel.copy()

        for i in range (self.wcsaxes):
            intermediate[i] -= self.crpix[i]

        intermediate = N.matrix (intermediate)
        intermediate = N.matrix (self.pc) * intermediate

        return intermediate.A

    def _lin_i2p (self, intermediate):
        """Do the linear projection to pixels.

        @param intermediate: one or more vectors (in pixels) in the form
            [[x]   or   [[x0  x1  x2]
             [y]]        [y0  y1  y2]]
        @type intermediate: ndarray (rank 2)

        @return: pixel coordinates:
            output = PC_matrix_inverse * intermediate + crpix
        @rtype: ndarray (rank 2)
        """

        intermediate = N.matrix (intermediate)
        pixel = N.matrix (self.inv_pc) * intermediate

        for i in range (self.wcsaxes):
            pixel[i] += self.crpix[i]

        return pixel.A

class celwcs (WCSbase):
    """This is a WCS object for celestial coordinates.

    This is based on the following two papers:
    Greisen, E. W., & Calabretta, M. R. 2002, A & A, 395, 1061 (Paper I),
    Calabretta, M. R. & Greisen, E. W., 2002, A & A, 395, 1077 (Paper II)

    The public methods are:
        wcs.info (verbose=False)
        world = wcs.frompixel (pixel)
        pixel = wcs.topixel (world, trim=True)
        i0 = wcs.longitudeAxis()
        i1 = wcs.latitudeAxis()
        keywords = wcs.getKeywords()
        ctype = wcs.getCtype()
        (crval, crpix, cd) = wcs.getRADec_WCS()
        wcs.setRADec_WCS (crval=None, crpix=None, cd=None)
        intermediate = wcs.pixelIntermediate (pixel)

    info() prints the values of the WCS parameters.
    frompixel() converts from pixel coordinates to world coordinates; the
        type and order of world coordinates is defined by the CTYPE keywords
        (see getCtype).
    topixel() converts from world coordinates to pixel coordinates.
    longitudeAxis() and latitudeAxis() return the axis number (zero indexed)
        for right ascension (or longitude) and declination (or latitude)
    getCtype() returns the array of CTYPE strings; this is needed in order
        to interpret the output of frompixel and to give the correct input
        to topixel.
    getRADec_WCS() returns some of the in-memory WCS keyword values, but
        only the elements relevant to right ascension and declination,
        and in that order.
    setRADec_WCS() overwrites some of the in-memory WCS keyword values;
        as with getRADec_WCS, the arrays must be 2-D, for RA and Dec.
    pixelIntermediate() converts from pixel coordinates to intermediate
        world coordinates, but the returned value includes only the two
        values xi & eta (tangent-plane coordinates), in that order.

    naxis is the dimension of a pixel coordinate vector, while wcsaxes is
    the dimension of a world coordinate vector.  wcsaxes must be at least
    2 if the WCS transformation includes celestial coordinates.  The
    default for naxis is 2, and the default for wcsaxes is naxis; either
    or both may be specified by header keywords (or dictionary keys).

    pixel is an array of length naxis (or it may optionally be padded
    with zeros to length wcsaxes.)  It may be rank one, in which case
    it is interpreted as a single coordinate vector, or it may be rank
    two, in which case it is expected to have shape (naxis, N), where
    N is the number of coordinate vectors.

    world is an array of length wcsaxes.  As with pixel, world may be
    either rank one or rank two.  Two elements of the world coordinate
    vector may be right ascension and declination (or other longitude/
    latitude pair).  If wcsaxes is greater than naxis, the other vector
    elements are linear functions of the pixel coordinates; one example
    would be wavelength for a long slit image.

    xi and eta are intermediate world coordinates for spherical coordinates
    (right ascension and declination, or longitude and latitude).  Like
    right ascension and declination, xi increases eastward and eta increases
    northward.

>>> wcs = {"wcsaxesz": 2, "mjd-obs": 999., "ctype1z": "RA", "ctype2z": "DEC", "crpix1z": 24., "crpix2z": 41., "crval1z": 75.5, "crval2z": 33.5, "cdelt1z": -0.02, "cdelt2z": 0.02, "pc1_1z": 0.5, "pc1_2z": 0.866025403784, "pc2_1z": -0.866025403784, "pc2_2z": 0.5}
>>> w = celwcs (wcs, alternate="z")
>>> pixel = N.zeros ((2,5), dtype=N.float64)
>>> pixel[0,0] = 23; pixel[1,0] = 20
>>> pixel[0,1] = 23; pixel[1,1] = 40
>>> pixel[0,2] = 23; pixel[1,2] = 60
>>> pixel[0,3] = 13; pixel[1,3] = 40
>>> pixel[0,4] = 33; pixel[1,4] = 40
>>> world = w.frompixel (pixel)
>>> print world
[[ 75.91445198  75.5         75.08362854  75.62016075  75.38031915]
 [ 33.29931296  33.5         33.69930082  33.67314641  33.32673807]]
>>> print w.topixel (world)
[[ 23.  23.  23.  13.  33.]
 [ 20.  40.  60.  40.  40.]]
>>> wcs = {"wcsaxes": 3, "ctype1": "WAVE", "ctype2": "RA---TAN", "ctype3": "DEC--TAN", "crpix1": 24., "crpix2": 41., "crval1": 5500., "crval2": 75.5, "crval3": 33.5, "cunit1": "Angstrom", "cdelt1": 0.4, "cdelt2": -0.02, "cdelt3": 0.02, "pc1_1": 1., "pc1_2": 0., "pc2_1": 0., "pc2_2": 0.866025403784, "pc3_1": 0., "pc3_2": 0.5}
>>> w = celwcs (wcs)
>>> world = w.frompixel (pixel)
>>> print world
[[ 5500.          5500.          5500.          5496.          5504.        ]
 [   75.91445198    75.5           75.08362854    75.5           75.5       ]
 [   33.29931296    33.5           33.69930082    33.5           33.5       ]]
>>> print w.topixel (world)
[[ 23.  23.  23.  13.  33.]
 [ 20.  40.  60.  40.  40.]]
>>> wcs = {"ctype1": "SOLY", "ctype2": "SOLX", "crpix1": 541.841247558594, "crpix2": 481.76220703125, "cunit1": "solRad", "cunit2": "solRad", "cd1_1": -0.0024600299075246, "cd1_2": 0.00107712228782475, "cd2_1": -0.0010762000456452, "cd2_2": -0.002455944661051}
>>> w = celwcs (wcs)
Warning:  longitude and latitude axes not found;
linear projection will be assumed.
>>> print w.frompixel ((540.841247558594, 480.76220703125))
[ 0.  0.]
>>> print w.topixel ((0., 0.))
[ 540.84124756  480.76220703]
>>> print w.frompixel ((541.841247558594, 481.76220703125))
[-0.00138291 -0.00353214]
>>> print w.topixel ((-0.0013829076197, -0.0035321447067))
[ 541.84124756  481.76220703]
>>> # ra and dec not in the usual order
>>> wcs = {"wcsaxes": 2, "radesys": "FK5", "ctype1": "DEC--TAN", "ctype2": "RA---TAN", "crpix1": 1094.5, "crpix2": 599.5, "crval1": 64.068867233514, "crval2": 258.033013388337, "cd1_1": 11., "cd1_2": 12., "cd2_1": 21., "cd2_2": 22.}
>>> w = celwcs (wcs)
>>> # the order of the results here is (dec, ra)
>>> print w.frompixel ((1103.5, 598.5))
[ 24.11057816   7.41302373]
>>> print w.frompixel ((1093.5, 608.5))
[ 23.86443351   8.67201834]
>>> # these methods return parameters in the order:  ra, dec
>>> (crval, crpix, cd) = w.getRADec_WCS()
>>> print crval
[ 258.03301339   64.06886723]
>>> print crpix
[ 1093.5   598.5]
>>> print cd
[[ 21.  22.]
 [ 11.  12.]]
>>> print w.pixelIntermediate ((1103.5, 598.5))
[ 210.  110.]
>>> # image array
>>> wcs = {"wcsnamez": "celestial", "ctype1z": "RA---TAN", "ctype2z": "DEC--TAN", "crpix1z": 257., "crpix2z": 513., "pc1_1z": 0.00011, "pv1_13z": 13.}
>>> w = celwcs (wcs, alternate="z")
>>> # bintable array
>>> wcs = {"wcax3z": 2, "wcsn3z": "celestial", "1cty3z": "RA---TAN", "2cty3z": "DEC--TAN", "equi3z": 1950., "1crp3z": 257., "2crp3z": 513., "11pc3z": 0.00011, "1v3_19z": 19.}
>>> w = celwcs (wcs, alternate="z", column=2)
>>> # pixel list (shorter keyword names)
>>> wcs = {"wcsn3z": "celestial", "tcty3z": "RA---TAN", "tcty4z": "DEC--TAN", "equi3z": 1950., "tcrv3z": 342.6, "tcrv4z": -75.3, "tcun3z": "deg", "tcun4z": "deg", "tcrp3z": 257., "tcrp4z": 513., "tcde3z": 0.00031, "tcde4z": 0.00041, "tp3_3z": 1.003, "tp4_4z": 1.004, "tv3_27z": 27., "tv4_43z": 43.,}
>>> w = celwcs (wcs, alternate="z", column=(2,3))
>>> # pixel list (longer keyword names)
>>> wcs = {"wcsn3z": "celestial", "tctyp3z": "RA---TAN", "tctyp4z": "DEC--TAN", "equi3z": 1950., "tcrvl3z": 342.6, "tcrvl4z": -75.3, "tcuni3z": "deg", "tcuni4z": "deg", "tcrpx3z": 257., "tcrpx4z": 513., "tcdlt3z": 0.00031, "tcdlt4z": 0.00041, "tpc3_3z": 1.003, "tpc4_4z": 1.004, "tpv3_27z": 27., "tpv4_43z": 43.,}
>>> w = celwcs (wcs, alternate="z", column=(2,3))
"""

    def __init__ (self, hdr, alternate=None, column=None):
        """constructor

        @param hdr: a set of FITS keywords and values (values are case
            sensitive, but keyword names are not)
        @type hdr: pyfits.Header object or dictionary
        @param alternate: a single letter, specifying an alternate WCS
            (default is blank, for primary WCS)
        @type alternate: string (one letter) or None
        @param column: for an image array, this should be None;
            for an array in a table column, this should be a single column
            number (zero indexed);
            for a tabulated list of pixels in two columns, this should be a
            list giving the zero-indexed column numbers for ctype1 and ctype2
            (nominally right ascension and declination)
        @type column: None, int, or list of two ints
        """

        # extract keywords using base class constructor
        WCSbase.__init__(self, hdr, alternate, column)

        # convert from degrees to radians, set attributes derived from keywords
        self._wcsSet()

    ### Beginning of public attributes.

    def info (self, verbose=False):
        """print info about the WCS

        @param verbose: if true, also print some of the more obscure info
        @type verbose: boolean
        """

        # This array is for converting RA and Dec from radians to degrees,
        # without affecting values for other axes (if any).
        RadiansToDegrees = N.ones (self.wcsaxes, dtype=N.float64)
        if self.longitude_axis is not None:
            RadiansToDegrees[self.longitude_axis] /= DEGREEStoRADIANS
        if self.latitude_axis is not None:
            RadiansToDegrees[self.latitude_axis] /= DEGREEStoRADIANS

        if self.alternate != "":
            print "alternate WCS '%s'" % self.alternate
        print "naxis = %d, wcsaxes = %d" % (self.naxis, self.wcsaxes)
        print "wcsname =", repr (self.wcsname)
        print "ctype =", self.ctype
        print "cunit =", self.cunit
        print "crval =", self.crval * RadiansToDegrees
        print "crpix =", self.crpix, " (zero indexed)"
        print "cdelt =", self.cdelt * RadiansToDegrees
        print "pc ="
        print self.pc
        if self.got_cd:
            print "cd ="
            print self.cd               # was not converted to radians
        if self.got_crota:
            print "crota =", self.crota * RadiansToDegrees
        if self.lonpole is not None:
            print "lonpole =", self.lonpole / DEGREEStoRADIANS,
            print " (native longitude of the celestial pole)"
        elif verbose:
            print "lonpole =", self.lonpole
        if self.latpole is not None:
            print "latpole =", self.latpole / DEGREEStoRADIANS,
            print " (native latitude of the celestial pole)"
        elif verbose:
            print "latpole =", self.latpole
        if self.theta_a is not None:
            print "theta_a =", self.theta_a / DEGREEStoRADIANS,
            print " (from pv keyword, for conic projection)"
            print "eta =", self.eta / DEGREEStoRADIANS,
            print " (from pv keyword, or default of 0)"
        if verbose:
            print "long_0 =", self.long_0/ DEGREEStoRADIANS,
            print " (native longitude of the fiducial point)"
            print "lat_0 =", self.lat_0/ DEGREEStoRADIANS,
            print " (native latitude of the fiducial point)"
            print "ra_0 =", self.ra_0 / DEGREEStoRADIANS,
            print " (celestial longitude of the fiducial point)"
            print "dec_0 =", self.dec_0 / DEGREEStoRADIANS,
            print " (celestial latitude of the fiducial point)"
            print "ra_p =", self.ra_p / DEGREEStoRADIANS,
            print " (celestial longitude of the native pole)"
            print "dec_p =", self.dec_p / DEGREEStoRADIANS,
            print " (celestial latitude of the native pole)"
        if self.mjd_obs is not None:
            print "mjd-obs =", self.mjd_obs
        print "equinox =", self.equinox
        print "radesys =", repr (self.radesys)
        for (i, m, value) in self.pv:
            print "pv%d_%d =" % (i+1, m), value
        for (i, m, value) in self.ps:
            print "ps%d_%d =" % (i+1, m), value
        if self.longitude_axis is not None:
            print "longitude axis =", self.longitude_axis, " (zero indexed)"
        if self.latitude_axis is not None:
            print "latitude axis =", self.latitude_axis, " (zero indexed)"
        print "map projection =", self.projection

    def frompixel (self, pixel):
        """Compute world coordinates from pixel coordinates.

        @param pixel: an array of one or more pixel positions
        @type pixel: sequence object

        @return: pixel vector(s) converted to world coordinates
        @rtype: ndarray

        If pixel contains just a single position vector, it can be written
        either as numpy.array ([x, y]) or as numpy.array ([[x], [y]]), i.e.
          [x  y] or [[x]
                     [y]]
        The function value will be in the same format (but will be padded
        with zeros if the vector dimension is less than WCSAXES).
        If pixel contains more than one position vector, it must be in
        the form:
          [[x0  x1  x2  x3]
           [y0  y1  y2  y3]]
        """

        pixel = N.array (pixel)
        if len (pixel.shape) == 1:
            single_vector = True        # pixel is a single vector
            pixel = N.reshape (pixel, (pixel.shape[0], 1))
        else:
            single_vector = False       # pixel is an array of positions

        # Pad the pixel coordinates with zeros, if necessary, to ensure
        # that the length of the coordinate vector is wcsaxes.
        pixel = self._padPixel (pixel)

        # A distortion correction (CPDIS) could be included at this point.

        # Do the linear part of the transformation.
        intermediate = self._lin_p2i (pixel)

        # A distortion correction (CQDIS) could be included at this point.

        # Convert from pixels to physical units.
        world = intermediate.copy()
        for i in range (self.wcsaxes):
            world[i] *= self.cdelt[i]
        del intermediate

        if self.projection != "linear":

            x = world[self.longitude_axis]
            y = world[self.latitude_axis]

            # include spherical projection
            (phi, theta) = self._inverse_fcn (x, y)

            (ra, dec) = self._nativeToCel (phi, theta)

            world[self.longitude_axis] = ra / DEGREEStoRADIANS
            world[self.latitude_axis] = dec / DEGREEStoRADIANS

        for i in range (self.wcsaxes):
            if i != self.longitude_axis and i != self.latitude_axis:
                world[i] += self.crval[i]

        if single_vector:
            world = N.reshape (world, (world.shape[0],))

        return world

    def topixel (self, world, trim=True):
        """Compute pixel coordinates from world coordinates.

        @param world: an array of one or more position vectors
        @type world: sequence object
        @param trim: if true and WCSAXES > NAXIS, trim zero padding from
            the output pixel coordinates before returning them
        @type trim: boolean

        @return: 'world' converted to pixels
        @rtype: ndarray

        If 'world' contains just a single position vector, it can be written
        either as numpy.array ([ra, dec]) or as numpy.array ([[ra], [dec]]),
          i.e. [ra  dec] or [[ra]
                             [dec]]
        The function value will be in the same format (but may be padded
        with zeros if trim is False).
        If 'world' contains more than one position vector, it must be in
        the form:
          [[ra0   ra1   ra2   ra3]
           [dec0  dec1  dec2  dec3]]
        """

        if self.inv_pc is None:
            raise ValueError, \
                  "Can't compute pixel coordinates because matrix is singular."

        intermediate = N.array (world)
        if len (intermediate.shape) == 1:
            single_vector = True        # world is a single vector
            intermediate = N.reshape (intermediate, (intermediate.shape[0], 1))
        else:
            single_vector = False       # world is an array of positions

        for i in range (self.wcsaxes):
            if i != self.longitude_axis and i != self.latitude_axis:
                intermediate[i] -= self.crval[i]

        if self.projection != "linear":

            ra = intermediate[self.longitude_axis] * DEGREEStoRADIANS
            dec = intermediate[self.latitude_axis] * DEGREEStoRADIANS

            (phi, theta) = self._celToNative (ra, dec)

            # include spherical projection
            (x, y) = self._projection_fcn (phi, theta)

            intermediate[self.longitude_axis] = x
            intermediate[self.latitude_axis] = y

        # Convert from physical units to pixels.
        for i in range (self.wcsaxes):
            intermediate[i] /= self.cdelt[i]

        # A distortion correction (CQDIS) could be included at this point.

        pixel = self._lin_i2p (intermediate)

        # A distortion correction (CPDIS) could be included at this point.

        if single_vector:
            pixel = N.reshape (pixel, (pixel.shape[0],))

        if trim:
            return pixel[0:self.naxis]          # the actual pixel coordinates
        else:
            return pixel

    def longitudeAxis (self):
        """Return the index for right ascension.

        @return: axis number (zero indexed) for right ascension (or longitude)
        @rtype: integer
        """

        return self.longitude_axis

    def latitudeAxis (self):
        """Return the index for right ascension.

        @return: axis number (zero indexed) for declination (or latitude)
        @rtype: integer
        """

        return self.latitude_axis

    def getCtype (self):
        """Return the CTYPE array.

        @return: CTYPE
        @rtype: array of character strings
        """

        return self.ctype

    def getKeywords (self):
        """Return the list of WCS keywords in hdr.

        @return: WCS keywords
        @rtype: list
        """

        keywords = copy.copy (self.keywords)
        keywords.sort()

        return keywords

    def getRADec_WCS (self):
        """Extract WCS parameters in (RA,Dec) order.

        @return: (crval, crpix, cd), in degrees, pixels, degrees per pixel.
            crval is the right ascension and declination at the reference pixel
            crpix is the reference pixel (one indexed)
            cd is the 2x2 CD matrix
        @rtype: tuple containing three two-element (or 2x2) ndarrays
        """

        if self.longitude_axis is None or self.latitude_axis is None:
            raise ValueError, "WCS does not include RA and Dec"

        crval = N.zeros (2, dtype=N.float64)
        crpix = N.zeros (2, dtype=N.float64)
        cd = N.identity (2, dtype=N.float64)

        crval[0] = self.crval[self.longitude_axis] / DEGREEStoRADIANS
        crval[1] = self.crval[self.latitude_axis]  / DEGREEStoRADIANS
        crpix[0] = self.crpix[0]
        crpix[1] = self.crpix[1]

        # Extract the elements pertaining to the spherical coordinates.
        i = self.i
        j = self.j
        cd[0,0] = self.cdelt[i] * self.pc[i,i] / DEGREEStoRADIANS
        cd[0,1] = self.cdelt[i] * self.pc[i,j] / DEGREEStoRADIANS
        cd[1,0] = self.cdelt[j] * self.pc[j,i] / DEGREEStoRADIANS
        cd[1,1] = self.cdelt[j] * self.pc[j,j] / DEGREEStoRADIANS
        if self.longitude_axis > self.latitude_axis:
            temp = cd[0].copy()
            cd[0] = cd[1].copy()
            cd[1] = temp.copy()

        return (crval, crpix, cd)

    def setRADec_WCS (self, crval=None, crpix=None, cd=None):
        """Set WCS parameters in (RA,Dec) order."""

        if self.longitude_axis is None or self.latitude_axis is None:
            raise ValueError, "WCS does not include RA and Dec"

        if crval is not None:
            crval = N.array (crval)
            if crval.shape[0] != 2:
                raise ValueError, "crval must be a two-element array"
            self.crval[self.longitude_axis] = crval[0] * DEGREEStoRADIANS
            self.crval[self.latitude_axis]  = crval[1] * DEGREEStoRADIANS

        if crpix is not None:
            crpix = N.array (crpix)
            if crpix.shape[0] != 2:
                raise ValueError, "crpix must be a two-element array"
            self.crpix[0] = crpix[0]
            self.crpix[1] = crpix[1]

        if cd is not None:
            cd = N.array (cd)
            if cd.shape != (2,2):
                raise ValueError, "cd must be a 2x2 array"
            # Compute the PC matrix and CDELT from CD.
            norm0 = cd[0,0]**2 + cd[0,1]**2
            norm1 = cd[1,0]**2 + cd[1,1]**2
            if norm0 <= 0. or norm1 <= 0:
                raise ValueError, "The CD matrix is singular."
            norm0 = math.sqrt (norm0)
            norm1 = math.sqrt (norm1)
            self.cdelt[self.longitude_axis] = norm0 * DEGREEStoRADIANS
            self.cdelt[self.latitude_axis]  = norm1 * DEGREEStoRADIANS
            tcd = cd.copy()
            if self.longitude_axis > self.latitude_axis:
                temp = tcd[0].copy()
                tcd[0] = tcd[1].copy()
                tcd[1] = temp.copy()
            self.pc[self.i,self.i] = tcd[0,0] / norm0
            self.pc[self.i,self.j] = tcd[0,1] / norm0
            self.pc[self.j,self.i] = tcd[1,0] / norm1
            self.pc[self.j,self.j] = tcd[1,1] / norm1

    def pixelIntermediate (self, pixel):
        """Compute intermediate world coordinates from pixel coordinates.

        @param pixel: an array of one or more pixel positions
        @type pixel: sequence object

        @return: pixel vector(s) converted to intermediate world coordinates
        @rtype: ndarray

        See frompixel() for a description of the input and output array
        formats.
        """

        if self.longitude_axis is None or self.latitude_axis is None:
            raise ValueError, "WCS does not include RA and Dec"

        pixel = N.array (pixel)
        if len (pixel.shape) == 1:
            single_vector = True        # pixel is a single vector
            pixel = N.reshape (pixel, (pixel.shape[0], 1))
        else:
            single_vector = False       # pixel is an array of positions

        # Pad the pixel coordinates with zeros, if necessary, to ensure
        # that the length of the coordinate vector is wcsaxes.
        pixel = self._padPixel (pixel)

        # A distortion correction (CPDIS) could be included at this point.

        # Do the linear part of the transformation.
        intermediate = self._lin_p2i (pixel)

        # A distortion correction (CQDIS) could be included at this point.

        # Convert from pixels to physical units.
        for i in range (self.wcsaxes):
            intermediate[i] *= self.cdelt[i]

        if single_vector:
            intermediate = N.reshape (intermediate, (intermediate.shape[0],))

        intermediate[self.longitude_axis] /= DEGREEStoRADIANS
        intermediate[self.latitude_axis] /= DEGREEStoRADIANS

        if single_vector:
            newshape = (2,)
        else:
            newshape = (2,intermediate.shape[1])
        temp = N.zeros (newshape, dtype=N.float64)
        temp[0] = intermediate[self.longitude_axis]
        temp[1] = intermediate[self.latitude_axis]

        return temp

    ### End of public attributes.

    def _initializeWcsAttributes (self):
        """Initialize WCS attributes to default values.

        We need to know the dimension of the WCS coordinate system before
        calling this method, because some of the attributes are arrays that
        depend on this size.
        """

        # self.naxis and self.wcsaxes were set by _findWcsAxes()
        self.ctype = []
        self.cunit = []
        for i in range (self.wcsaxes):
            self.ctype.append ("PIXEL")
            self.cunit.append ("")
        self.crval = N.zeros (self.wcsaxes, dtype=N.float64)
        self.crpix = N.zeros (self.wcsaxes, dtype=N.float64)
        self.cdelt = N.ones (self.wcsaxes, dtype=N.float64)
        self.pc = N.identity (self.wcsaxes, dtype=N.float64)
        self.inv_pc = None
        self.cd = N.zeros ((self.wcsaxes,self.wcsaxes), dtype=N.float64)
        self.crota = N.zeros (self.wcsaxes, dtype=N.float64)
        self.crder = N.zeros (self.wcsaxes, dtype=N.float64)
        self.csyer = N.zeros (self.wcsaxes, dtype=N.float64)
        # keywords lonpole and latpole are included below.
        self.wcsname = ""
        self.mjd_obs = None
        self.equinox = None
        self.radesys = None
        self.pv = []
        self.npv = 0
        self.ps = []
        self.nps = 0

        # These are boolean flags to indicate whether we have keywords
        # of this type.
        self.got_pc = False
        self.got_cd = False
        self.got_cdelt = False
        self.got_crota = False

        # These specify the axis numbers for RA and Dec (or other
        # longitude/latitude pair).
        self.longitude_axis = None
        self.latitude_axis = None
        self.i = None
        self.j = None

        # See Table 1 of Calabretta & Greisen (Paper II).
        # Celestial longitude and latitude of the fiducial point.
        self.ra_0 = None                        # CRVAL<longitude_axis>
        self.dec_0 = None                       # CRVAL<latitude_axis>
        # Celestial longitude and latitude of the native pole.
        self.ra_p = None
        self.dec_p = None
        # Native longitude and latitude of the fiducial point.
        self.long_0 = None                      # default, or PV1_1
        self.lat_0 = None                       # default, or PV1_2
        # Native longitude and latitude of the celestial pole (long_p, lat_p).
        self.lonpole = None                     # LONPOLE, or PV1_3
        self.latpole = None                     # LATPOLE, or PV1_4
        # These are used for conic projections.
        self.theta_a = None                     # PV2_1 (latitude)
        self.eta = 0.                           # PV2_2

        # This is the spherical map projection (e.g. "TAN").
        self.projection = "linear"

    def _wcsSet (self):
        """Do some computations and interpretations of the WCS info."""

        if self.got_cdelt:
            for i in range (self.wcsaxes):
                if self.cdelt[i] == 0.:
                    raise ValueError, "CDELT must not be zero."

        if self.got_cd:
            if self.got_pc:
                print "Warning:  Both PC and CD keywords were found;"
                print "the CD keywords will be ignored, "\
                      "and PC will be used instead."
                self.got_cd = False
            elif self.got_cdelt:
                print "Warning:  Both CDELT and CD keywords were found;"
                print "the CDELT keywords will be recomputed from CD."

        if self.got_crota:
            if self.got_pc:
                key = "PC"
            elif self.got_cd:
                key = "CD"
            if self.got_pc or self.got_cd:
                print "Warning:  Both CROTA and", key, "keywords were found;"
                print "the CROTA keywords will be ignored."
                self.got_crota = False

        if self.lonpole is not None:
            self.lonpole *= DEGREEStoRADIANS
        if self.latpole is not None:
            self.latpole *= DEGREEStoRADIANS

        r1 = re.compile (r"RA|GLON|ELON|HLON|SLON.*")
        r2 = re.compile (r"DEC|GLAT|ELAT|HLAT|SLAT.*")
        for i in range (len (self.ctype)):
            if r1.match (self.ctype[i]) is not None:
                self.longitude_axis = i
            elif r2.match (self.ctype[i]) is not None:
                self.latitude_axis = i
        messages = []
        if self.longitude_axis is None:
            messages.append ("longitude")
        if self.latitude_axis is None:
            messages.append ("latitude")
        if len (messages) > 0:
            if len (messages) == 1:
                print "Warning: ", messages[0], "axis not found;"
            if len (messages) == 2:
                print "Warning:  longitude and latitude axes not found;"
            print "linear projection will be assumed."
            self.projection = "linear"

        if self.longitude_axis is not None and self.latitude_axis is not None:
            self.i = min (self.longitude_axis, self.latitude_axis)
            self.j = max (self.longitude_axis, self.latitude_axis)

        # Convert ra & dec parameters from degrees to radians.
        if self.longitude_axis is not None:
            i_ra = self.longitude_axis
            self.crval[i_ra] *= DEGREEStoRADIANS
            self.ra_0 = self.crval[i_ra]
            self.cdelt[i_ra] *= DEGREEStoRADIANS
            self.crder[i_ra] *= DEGREEStoRADIANS
            self.csyer[i_ra] *= DEGREEStoRADIANS
            self.crota[i_ra] *= DEGREEStoRADIANS
        if self.latitude_axis is not None:
            i_dec = self.latitude_axis
            self.crval[i_dec] *= DEGREEStoRADIANS
            self.dec_0 = self.crval[i_dec]
            self.cdelt[i_dec] *= DEGREEStoRADIANS
            self.crder[i_dec] *= DEGREEStoRADIANS
            self.csyer[i_dec] *= DEGREEStoRADIANS
            self.crota[i_dec] *= DEGREEStoRADIANS

        if self.got_cd:
            # Compute the PC matrix and CDELT from CD.
            for i in range (self.cd.shape[0]):
                norm2 = 0.
                for j in range (self.cd.shape[1]):
                    norm2 += (self.cd[i,j]**2)
                if norm2 <= 0.:
                    raise ValueError, "The CD matrix is singular."
                norm = math.sqrt (norm2)
                self.cdelt[i] = norm
                if i == self.longitude_axis or i == self.latitude_axis:
                    self.cdelt[i] *= DEGREEStoRADIANS
                self.pc[i] = self.cd[i] / norm
        elif self.got_crota and self.longitude_axis is not None and \
                                self.latitude_axis is not None:
            # Compute the PC matrix from CROTA and CDELT.
            # See section 6.1 of Calabretta & Greisen (Paper II).
            i = self.i
            j = self.j
            crota = self.crota[j]
            ratio = self.cdelt[j] / self.cdelt[i]
            self.pc[i,i] =  math.cos (crota)
            self.pc[i,j] = -math.sin (crota) * ratio
            self.pc[j,i] =  math.sin (crota) / ratio
            self.pc[j,j] =  math.cos (crota)

        try:
            # compute the inverse of the PC matrix
            pc = N.matrix (self.pc)
            self.inv_pc = pc.I.A
        except:
            print "Warning:  the linear transformation matrix is singular."
            self.inv_pc = None

        # The defaults for radesys and equinox depends on which keywords
        # were found and what their values are, so we can't assign defaults
        # until we've read all the keywords.
        if self.radesys is None:
            if self.equinox is None:
                self.equinox = 2000.
                self.radesys = "ICRS"
            elif self.equinox < 1984.:
                self.radesys = "FK4"
            else:
                self.radesys = "FK5"
        elif self.equinox is None:
            if self.radesys[0:3] == "FK4":
                self.equinox = 1950.
            else:
                self.equinox = 2000.

        if self.longitude_axis is not None:
            self._setSphProj()          # spherical projection attributes

        if self.longitude_axis is not None and self.latitude_axis is not None:
            self._setFromPV()           # use PV keywords, if specified
            self._defaultLongLat()      # default lonpole, latpole, etc
            self._setRaDec_p()          # set ra_p, dec_p attributes

    def _setSphProj (self):
        """Set the spherical projection name and function attributes."""

        ctype_i = self.ctype[self.longitude_axis]
        parts = ctype_i.split ("-")
        proj = parts[-1]
        if proj == "" or proj == parts[0]:
            self.projection = "TAN"         # default spherical projection
        elif proj in SPHERICAL_MAP_PROJ:
            self.projection = proj
        else:
            print "Warning:  '%s' projection is not recognized;" % (proj)
            print "tangent projection will be assumed."
            self.projection = "TAN"
        # Set the attributes for the spherical projection functions.
        p_attrib = self.projection + "_projection"
        invp_attrib = self.projection + "_inverse"
        try:
            self._projection_fcn = self.__getattribute__ (p_attrib)
            self._inverse_fcn = self.__getattribute__ (invp_attrib)
        except AttributeError:
            raise RuntimeError, \
                  self.projection + " projection is not implemented"

    def _setFromPV (self):
        """Set parameters from the PV keywords.

        Note that lonpole or latpole set by a keyword of that name may be
        overridden by a PV keyword.
        """

        for k in range (self.npv):
            (i, m, value) = self.pv[k]
            if i == self.longitude_axis:
                # set long_0, lat_0, LONPOLE or LATPOLE
                if m == 1:
                    self.long_0 = value * DEGREEStoRADIANS
                elif m == 2:
                    self.lat_0  = value * DEGREEStoRADIANS
                elif m == 3:
                    self.lonpole = value * DEGREEStoRADIANS
                elif m == 4:
                    self.latpole = value * DEGREEStoRADIANS
            elif i == self.latitude_axis:
                # for conic projections
                if m == 1:
                    self.theta_a = value * DEGREEStoRADIANS
                    if value == 0.:
                        print "Warning:  PV%d_%d = 0 will cause problems " \
                              "for conic projections" % (i+1, m)
                elif m == 2:
                    self.eta  = value * DEGREEStoRADIANS

    def _defaultLongLat (self):
        """Assign defaults for long_0, lat_0, lonpole.

        Note that no default is assigned for latpole.
        """

        if self.long_0 is None:
            self.long_0 = 0.

        if self.lat_0 is None:
            if self.projection in ZENITHAL_PROJ:
                self.lat_0 = math.pi / 2.
            elif self.projection in CYLINDRICAL_PROJ:
                self.lat_0 = 0.
            elif self.projection in CONIC_PROJ:
                if self.theta_a is None:
                    raise ValueError, "theta_a is undefined (use a PV keyword)"
                self.lat_0 = self.theta_a
            else:
                self.lat_0 = 0.

        if self.lonpole is None:
            if self.dec_0 >= self.lat_0:
                self.lonpole = 0.
            else:
                self.lonpole = math.pi

    def _setRaDec_p (self):
        """Assign values for the celestial coordinates of the native pole."""

        if self.projection in ZENITHAL_PROJ:
            self.ra_p = self.ra_0
            self.dec_p = self.dec_0
        else:
            # Assign a value for latpole in case we need it.
            latpole = self.latpole
            if latpole is None:
                if self.projection in ZENITHAL_PROJ:
                    latpole = self.dec_0
                elif self.projection in CYLINDRICAL_PROJ:
                    latpole = math.pi / 2.
                elif self.projection in CONIC_PROJ:
                    latpole = math.pi / 2.
                else:
                    latpole = self.dec_0
            delta_long = self.lonpole - self.long_0
            sqrt_arg = 1. - math.cos (self.lat_0)**2 * math.sin (delta_long)**2
            if sqrt_arg <= TINY:
                self.dec_p = latpole
            else:
                temp1 = math.atan2 (math.sin (self.lat_0),
                            math.cos (self.lat_0) * math.cos (delta_long))
                temp2 = math.acos (math.sin (self.dec_0) / math.sqrt (sqrt_arg))
                dec1 = temp1 + temp2
                dec2 = temp1 - temp2
                if dec1 >= -math.pi / 2. and dec1 <= math.pi / 2. and \
                   dec2 >= -math.pi / 2. and dec2 <= math.pi / 2.:
                    if abs (dec1 - latpole) <= abs (dec2 - latpole):
                        self.dec_p = dec1
                    else:
                        self.dec_p = dec2
                elif dec1 >= -math.pi / 2. and dec1 <= math.pi / 2.:
                    self.dec_p = dec1
                elif dec2 >= -math.pi / 2. and dec2 <= math.pi / 2.:
                    self.dec_p = dec2
                else:
                    print "Warning:  can't compute declination " \
                          "at the native pole (90 degrees used)"
                    self.dec_p = math.pi / 2.

            if math.cos (self.dec_p) < TINY:
                if self.dec_p > 0:
                    self.ra_p = self.ra_0 + delta_long - math.pi
                else:
                    self.ra_p = self.ra_0 - delta_long
            else:
                sin_da = math.sin (delta_long) * math.cos (self.lat_0) / \
                         math.cos (self.dec_0)
                cos_da = (math.sin (self.lat_0) - \
                          math.sin (self.dec_p) * math.sin (self.dec_0)) / \
                         (math.cos (self.dec_p) * math.cos (self.dec_0))
                da = math.atan2 (sin_da, cos_da)
                self.ra_p = self.ra_0 - da

    def _nativeToCel (self, phi, theta):
        """Convert from native spherical coordinates to celestial."""

        sin_theta = N.sin (theta)
        cos_theta = N.cos (theta)
        sin_dphi  = N.sin (phi - self.lonpole)
        cos_dphi  = N.cos (phi - self.lonpole)
        sin_dec_p = N.sin (self.dec_p)
        cos_dec_p = N.cos (self.dec_p)

        x_temp = sin_theta * cos_dec_p - cos_theta * sin_dec_p * cos_dphi
        y_temp = -cos_theta * sin_dphi
        ra = N.arctan2 (y_temp, x_temp) + self.ra_p
        del x_temp, y_temp
        twopi = 2. * N.pi
        ra = N.remainder (ra+twopi, twopi)

        z_temp = sin_theta * sin_dec_p + cos_theta * cos_dec_p * cos_dphi
        dec = N.arcsin (z_temp)

        return (ra, dec)

    def _celToNative (self, ra, dec):
        """Convert from celestial coordinates to native spherical coords."""

        dra = ra - self.ra_p

        sin_dra  = N.sin (dra)
        cos_dra  = N.cos (dra)
        sin_dec = N.sin (dec)
        cos_dec = N.cos (dec)
        sin_dec_p = N.sin (self.dec_p)
        cos_dec_p = N.cos (self.dec_p)

        x_temp = sin_dec * cos_dec_p - cos_dec * sin_dec_p * cos_dra
        y_temp = -cos_dec * sin_dra
        phi = N.arctan2 (y_temp, x_temp) + self.lonpole
        del x_temp, y_temp

        z_temp = sin_dec * sin_dec_p + cos_dec * cos_dec_p * cos_dra
        theta = N.arcsin (z_temp)

        return (phi, theta)

def _test():
    import doctest
    import celwcs
    return doctest.testmod (celwcs)

if __name__ == "__main__":
    _test()
