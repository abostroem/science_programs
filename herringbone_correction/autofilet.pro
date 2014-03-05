PRO AUTOFILET, image, result, fpspec, EXTEN=exten,$
		XORDER=xorder, YORDER=yorder, FITX=fitx, SCUT=scut, SCNB=scnb,$
		NREJ=nrej, NGROW=ngrow, FREQ1=freq1, FREQ2=freq2,$
		TEXTLOG=textlog, DEBUG=debug

;+ ---------------------------------------------------------------------
; Program: AUTOFILET v1.1    Author: Rolf A. Jansen (ASU) -- Aug 30 2007
;
; Procedure to fit the herring-bone pattern in the read-noise of the HST
; STIS side-2 electronics CCD frames (or similar pick-up noise in CCD
; data from other telescopes), by optimized Fourier filtering of the
; frequencies responsible for producing the pattern.  The pattern is fit
; in the Fourier transformed time-series constructed from a residuals
; image, instead of from the science image itself.  This approach avoids
; ringing due to genuine high-frequency non-periodic signals (objects,
; spectral lines, cosmic ray events).  Output is a 2-D FITS format image
; containing only the herring-bone pattern and a 2-D pattern-subtracted
; science image. To remove the herring-bone pattern from a science image
; that is part of a multi-layer image set [SCI,ERR,DQ], simply subtract
; the herring-bone image.
;
; This procedure builds on 'stisnoise.pro' by Thomas M. Brown (ISR STIS
; 2001-005) and revision thereof by Nick Collins, Rolf Jansen & Rogier
; Windhorst (GSFC, dd. Jan 11, 2002).
;
; Dr. R.A. Jansen                            e-mail: Rolf.Jansen@asu.edu
; Arizona State University                   phone:  (480) 727-7119
; P.O. Box 871404 / Tyler Mall PSF-686
; TEMPE, AZ 85287-1404, USA
;
; ----------------------------------------------------------------------
; Input:   image   (string)   - CCD image frame (STIS: IMAGE,SCI frame)
; Output:  result  (string)   - Cleaned (pattern subtracted) input image
;          fpspec  (string)   - Output frequency power spectrum table
;
; Optional input (and their defaults):
;          exten   (int)   (1) - FITS extension containing the SCI image
;          xorder  (int)  (15) - Order of spline3 fit along rows
;          yorder  (int)   (5) - Order of spline3 fit along columns
;          fitx    (int)   (0) - First fit along rows instead of columns
;          scut    (real) (3.) - (Raw) sigma-cut for genuine signal
;          scnb    (real) (.5) - Sigma-cut for neighbor of rejected pixel
;          nrej    (int)   (3) - Number of rejection iterations
;          ngrow   (int)   (1) - Rejection growing radius in pixels
;          freq1   (real)(1.e4)- Lower limit of frequency interval [Hz]
;          freq2   (real)(2.e4)- Upper limit of frequency interval [Hz]
;          textlog (string)("autofilet.log") - Text log file
;          debug   (int)   (0) - Debugging flag for extra verbosity
;
; Example:
;          IDL> autofilet, "o6dc9b040_raw.fits","o6dc9b040_rwc.fits",$
;               "o6dc9b040_rwh.fits","o6dc9b040_powspec.dat",$
;               xorder=15,yorder=4,nrej=5
;
; ----------------------------------------------------------------------
; USES splfit.pro, splinefit.pro, readfits.pro, writefits.pro, sxpar.pro
;      get_date.pro, sxaddpar.pro, sxdelpar.pro, rfunct.pro
; ----------------------------------------------------------------------
;
; Revision History:
;          v0.1 -- May  1 2002 (RAJ) Extracted relevant routines from
;                  'stisnoise.pro' by T.M. Brown and from the revision
;                  thereof by Collins, Jansen & Windhorst (Jan 11, 2002)
;          v0.2 -- May  3 2002 (RAJ) Improved automated peak finding.
;          v0.3 -- May  6 2002 (RAJ) Data model fitting step implemented.
;          v0.4 -- May  8 2002 (RAJ) Replaced multiplicative suppression
;                  (windowing function) by substitution by white noise.
;          v0.5 -- May  9 2002 (RAJ) Tweaked prime number for FFT's
;          v0.9 -- May 13 2002 (RAJ) Fully functional stable version.
;          v1.0 -- Aug  8 2002 (RAJ) Updated default program parameters;
;                  output image headers now conform FITS standard, added
;                  keywords HPATCORR and HISTORY to the output headers.
;          v1.1 -- Aug 20-30 2007 (RAJ) IDL6.3; Removed xfit/yfit options,
;                  removed splinefit plot window and removed frequency/
;                  magnitude output arrays, but added output of frequency
;                  power spectrum table for use in HST-AR-11258 program.
;                  Added check for anomalous patterns. Added fitx option
;                  to allow fitting along rows instead of columns first.
;
;- ----------------------------------------------------------------------
; Permission to freely use, copy, modify, and distribute this software
; and its documentation without fee is hereby granted for non-military
; and non-commercial purposes, provided that the authors are acknowledged
; where due, the revision history is updated and the following disclaimer
; of warranty appears in all copies.
;
; DISCLAIMER:
;
; THIS PROGRAM IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND,
; EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
; TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY
; IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
; PURPOSE, AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE
; DOCUMENTATION WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE
; SOFTWARE WILL BE ERROR FREE.  IN NO EVENT SHALL THE AUTHOR(S) BE LIABLE
; FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT,
; SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR
; IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON
; WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS
; SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS
; WAS SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE
; SOFTWARE OR SERVICES PROVIDED HEREUNDER.
;
; ----------------------------------------------------------------------

progver='AUTOFILET v1.1 - R.A. Jansen - (IDL 6.3)'

; CHECK THE NUMBER OF COMMAND LINE ARGUMENTS...
IF (n_params(0) LT 3) THEN BEGIN
  print,'SYNTAX:  AUTOFILET, "image", "result", "fpspec", [exten=#],'
  print,'  where "image" is the name of a FITS file, "result" the name of the'
  print,'  cleaned output image, and "fpspec" the name of the frequency power'
  print,'  spectrum table (ASCII), and exten=# specifies the FITS extension'
  print,'  containing the science image (default=1).'
  print,'FITTING OPTIONS:'
  print,'  [xorder=#], [yorder=#], [/fitx]'
  print,'PIXEL REJECTION OPTIONS:'
  print,'  [scut=#.#], [scnb=#.#], [nrej=#], [ngrow=#]'
  print,'PATTERN FREQUENCY OPTIONS:'
  print,'  [freq1=#####.#], [freq2=#####.#]           (search interval [Hz])'
  print,'TEXT LOG FILE:'
  print,'  [textlog="aaaaaaaa"]'
  RETURN
ENDIF

; PROGRAM DEFAULTS...
IF (NOT keyword_set(exten))   THEN exten   = 1     ; FITS extension
IF (NOT keyword_set(xorder))  THEN xorder  = 15    ; order of fit in x
IF (NOT keyword_set(yorder))  THEN yorder  = 5     ; order of fit in y
IF (NOT keyword_set(fitx))    THEN fitx    = 0     ; fit along rows first
IF (NOT keyword_set(scut))    THEN scut    = 3.0   ; sigma clipping
IF (NOT keyword_set(scnb))    THEN scnb    = 0.5   ; neighbor rejected pixel
IF (NOT keyword_set(nrej))    THEN nrej    = 3     ; rejection iterations
IF (NOT keyword_set(ngrow))   THEN ngrow   = 1     ; growing radius
IF (NOT keyword_set(freq1))   THEN freq1   = 1.0e4 ; lower frequency limit
IF (NOT keyword_set(freq2))   THEN freq2   = 2.0e4 ; upper frequency limit
IF (NOT keyword_set(textlog)) THEN textlog ="autofilet.log"
IF (NOT keyword_set(debug))   THEN debug   = 0     ; debugging level (0|1)

; OTHER DEFAULTS...
dc = 0		; set power in first frequency bin (0 frequency) to zero

; CHECK INPUT PARAMETERS (assume negative values to be typos)...
IF (exten LT 0)  THEN exten = abs(exten)
IF (xorder LT 0) THEN xorder = abs(xorder)
IF (yorder LT 0) THEN yorder = abs(yorder)
IF (scut LT 0.0) THEN scut = abs(scut)
IF (scnb LT 0.0) THEN scnb = abs(scnb)
IF (nrej LT 0)   THEN nrej = abs(nrej)
IF (ngrow LT 0)  THEN ngrow = abs(ngrow)
IF (freq1 LT 0.) THEN freq1 = abs(freq1)
IF (freq2 LT 0.) THEN freq2 = abs(freq2)
IF (scnb GE scut) THEN BEGIN
   print,'WARNING: the sigma clipping criterion for neighbors of rejected pixels'
   print,'  should be more relaxed than that for finding significantly deviant'
   print,'  pixels. Recommended choice for "scnb" is ~0.10-0.30*"scut", depending'
   print,'  on the details of your data.'
   print,'  Default values are: scut = 3.0 sigma, scnb = 0.5 sigma.'
   RETURN
ENDIF


; ----------------------------------------------------------------------

T0  = SYSTIME(1)
T0a = SYSTIME(0)

; PADD FILE NAME EXTENSION, IF MISSING...
imtyp = ".fits"
cext  = STRMID(image,(STRPOS(image,'.',/REVERSE_SEARCH)+1))
; Valid image format extensions are <=4 characters in length; if more,
; assume it is part of the rootname and that the extension is missing:
IF (STRLEN(cext) LE 4 AND cext NE "fits") THEN BEGIN
   print,'ERROR: image extension ".'+cext+'" is not ".fits"! Exit...'
   RETURN
ENDIF
IF (STRLEN(cext) GT 4) THEN image = image+imtyp

; For output images, append ".fits" if not already the extension:
cext = STRMID(result,STRPOS(result,'.',/REVERSE_SEARCH)+1)
IF (cext NE "fits") THEN result = result+imtyp

; For the frequency power spectrum output table, append '.dat' if no obvious
; extension is provided:
cext = STRMID(fpspec,STRPOS(fpspec,'.',/REVERSE_SEARCH)+1)
IF (cext NE "dat" AND cext NE "txt") THEN fpspec = fpspec+'.dat'


; READ THE FITS PRIMARY HEADER UNIT, AND CHECK WHETHER "EXTEND == T" ...
phu  = READFITS(image,primehead,exten=0,/silent)
xtnd = STRTRIM(SXPAR(primehead,'EXTEND'))
; READ MJD AT THE START OF THE 1ST EXPOSURE FROM THE PRIMARY HEADER UNIT...
tmjd = STRTRIM(SXPAR(primehead,'TEXPSTRT'))

; READ THE (SPECIFIED EXTENSION IN THE) FITS FILE...
IF xtnd THEN BEGIN
   img  = READFITS(image,extenhead,exten=exten,/silent)
   xtn  = '['+STRTRIM(SXPAR(extenhead,'EXTNAME'))+']'
   cxtn = '['+STRTRIM(exten,1)+'="'+STRTRIM(SXPAR(extenhead,'EXTNAME'))+'"]'
ENDIF ELSE BEGIN
   img  = READFITS(image,primehead,/silent)
   xtn  = ''
ENDELSE

; RETRIEVE SOME INFORMATION FROM THE HEADER...
telescop = STRTRIM(SXPAR(primehead,'TELESCOP'))
instrume = STRTRIM(SXPAR(primehead,'INSTRUME'))

; DEFINE PHYSICAL CHARACTERISTICS OF THE CCD DEVICE...
IF (telescop EQ "HST" AND instrume EQ "STIS") THEN BEGIN
   ; STIS/CCD (Side-2 electronics) aboard HST.
   amp  = STRTRIM(SXPAR(primehead,'CCDAMP'))
   targ = STRTRIM(SXPAR(primehead,'TARGNAME'))
   gain = STRTRIM(SXPAR(primehead,'CCDGAIN'))
   ; Note: 19 leading + 19 trailing physical columns of serial overscan
   ; and 20 trailing virtual rows of parallel overscan. 
   pst = 640.0	 ; parallel shift time (microsec)
   sst = 22.0	 ; serial shift period (microsec)
   nc0 = 1062	 ; number of columns in raw data frame
   nr0 = 1044	 ; number of rows in raw data frame
   nxy = 1024	 ; number of columns and rows in flatfielded data
   nos = 19	 ; number of physical overscan columns
   pps = pst/sst ; number of serial shifts in a parallel shift interval
ENDIF
IF (telescop EQ "TILLINGHAST" AND instrume EQ "FAST") THEN BEGIN
   ; Original FAST CCD (FAST1) at the FLWO 60" telescope (1994-'96)
   ; Ask Steve Amato <samato@cfa.harvard.edu> for times pst and sst.
   amp  = ""
   targ = STRTRIM(SXPAR(primehead,'OBJECT'))
   gain = STRTRIM(SXPAR(primehead,'GAIN'))
   pst = XXX.X   ; of order 650 musec too?
   sst = XX.X    ; of order 20 musec too?
   nc0 = 2720	 ; number of columns in raw data frame
   nr0 = 512	 ; number of rows in raw data frame
   nxx = 2688    ; number of columns in trimmed data
   nyy = 480	 ; number of rows in trimmed data
   nos = 32	 ; number of physical overscan columns
   pps = pst/sst
ENDIF
; replicate and adapt the following IF-block to the specifics of data
; from other telescopes/instruments or ccd-controlers ...
IF (telescop EQ "AAAAAA" AND instrume EQ "BBBB") THEN BEGIN
   amp  = ""
   targ = STRTRIM(SXPAR(primehead,'OBJECT'))
   gain = STRTRIM(SXPAR(primehead,'GAIN'))
   pst = XXX.X   ; parallel shift time (microsec)
   sst = XX.X    ; serial shift period (microsec)
   nc0 = XXXX	 ; number of columns in raw data frame
   nr0 = XXXX	 ; number of rows in raw data frame
   nxx = XXXX    ; number of columns in trimmed data (if nxx != nyy)
   nyy = XXXX    ; number of rows in trimmed data (if nxx != nyy)
   nxy = XXXX    ; number of trimmed columns and rows (if nxx == nyy)
   nos = XX	 ; number of physical overscan columns
   pps = pst/sst
ENDIF

; DETERMINE THE IMAGE SIZE (number of columns and lines)...
image_sz = SIZE(img)
ncol = image_sz[1]
nrow = image_sz[2]

; OPEN TEXT LOGFILE...
OPENW,25,textlog,/APPEND

; PRINT OPENING MESSAGE...
printf,25,' '
printf,25,'# AUTOFILET: '+SYSTIME(/UTC)
printf,25,'  image=',image,cxtn,'[',ncol,',',nrow,']: ',targ,$
      '  amp="',amp,'"  gain=',gain,f='(a,a,a,a,i4,a,i4,a,a,a,a,a,f4.2)'
print,' '
print,'# AUTOFILET: '+SYSTIME(/UTC)
print,'  image=',image,cxtn,'[',ncol,',',nrow,']: ',targ,$
      '  amp="',amp,'"  gain=',gain,f='(a,a,a,a,i4,a,i4,a,a,a,a,a,f4.2)'
FLUSH, 25

; ----------------------------------------------------------------------

; CHECKS/PRE-PROCESSING SPECIFIC TO HST/STIS IMAGES...
IF (instrume EQ "STIS") THEN BEGIN
   ; only run this routine on science extensions...
   IF (xtn NE "[SCI]") THEN BEGIN
      printf,25,'ERROR: not a science [SCI] extension!'
      print,'ERROR: not a science [SCI] extension!'
      FLUSH, 25
      RETURN
   ENDIF

   ; the image size should be 1062x1044 or 1024x1024 pixels...
   IF (ncol EQ nc0 AND nrow EQ nr0) THEN BEGIN
      type = 'raw'
   ENDIF ELSE BEGIN
      IF (ncol EQ nxy AND nrow EQ nxy) THEN BEGIN
         type = 'flt'
      ENDIF ELSE BEGIN
         printf,25,'ERROR: neither 1062x1044 (raw) nor 1024x1024 (flt)!'
         print,'ERROR: neither 1062x1044 (raw) nor 1024x1024 (flt)!'
         FLUSH, 25
         RETURN
      ENDELSE
   ENDELSE

   ; padd trimmed data with fake 'overscan' strips...
   IF (type EQ "flt") THEN BEGIN
      printf,25,'  padding ',image,' with "overscan" strips...'
      print,'  padding ',image,' with "overscan" strips...'
      FLUSH, 25
      tmp = FLTARR(nc0,nxy)
      FOR row=0,nxy-1 DO tmp[*,row] = MEDIAN(img[*,row],/EVEN)
      tmp[nos:(nc0-nos-1),*] = img
      img = tmp
      ncol = nc0
   ENDIF

   ; rotate/transpose frame such that it is in read-out order...
   IF (amp EQ "A") THEN BEGIN
      ; data read thru pre-amp A -> leave as is:
   ENDIF
   IF (amp EQ "B") THEN BEGIN
      ; data read thru pre-amp B -> flip in x:
      ;print,'  flipping "',image,'" in x...'
      img = REVERSE(TEMPORARY(img),1)
   ENDIF
   IF (amp EQ "C") THEN BEGIN
      ; data read thru pre-amp C -> flip in y:
      ;print,'  flipping "',image,'" in y...'
      img = REVERSE(TEMPORARY(img),2)
   ENDIF
   IF (amp EQ "D") THEN BEGIN
      ; data read thru nominal pre-amp D -> rotate 180 degrees:
      ;print,'  rotating "',image,'" over 180 degrees...'
      img = ROTATE(TEMPORARY(img),2)
   ENDIF
   IF (amp NE "A" AND amp NE "B" AND amp NE "C" AND amp NE "D") THEN BEGIN
      ; if CCDAMP keyword is missing:
      printf,25,'ERROR: no amplifier specified in header!'
      print,'ERROR: no amplifier specified in header!'
      FLUSH, 25
      RETURN
   ENDIF
ENDIF

; CHECKS/PRE-PROCESSING SPECIFIC TO FLWO60"/FAST1 IMAGES...
IF (instrume EQ "FAST") THEN BEGIN
   ; the image size should be 2720x512 or 2688x480 pixels...
   IF (ncol EQ nc0 AND nrow EQ nr0) THEN BEGIN
      type = 'raw'
   ENDIF ELSE BEGIN
      IF (ncol EQ nxx AND nrow EQ nyy) THEN BEGIN
         type = 'trimmed'
      ENDIF ELSE BEGIN
         printf,25,'ERROR: neither 2720x512 (raw) nor 2688x480 (trimmed)! Perhaps'
         printf,25,'       the CCD was binned in y on read-out. If so, edit script.'
         print,'ERROR: neither 2720x512 (raw) nor 2688x480 (trimmed)! Perhaps'
         print,'       the CCD was binned in y on read-out. If so, edit script.'
         FLUSH, 25
         RETURN
      ENDELSE
   ENDELSE

   ; padd trimmed data with fake 'overscan' strips...
   IF (type EQ "trimmed") THEN BEGIN
      printf,25,'  padding ',image,' with "overscan" strips...'
      print,'  padding ',image,' with "overscan" strips...'
      FLUSH, 25
      tmp = FLTARR(nc0,nyy)
      FOR row=0,nyy-1 DO tmp[*,row] = MEDIAN(img[*,row],/EVEN)
      tmp[nos:(nc0-nos-1),*] = img
      img = tmp
      ncol = nc0
   ENDIF
ENDIF

; CHECKS/PRE-PROCESSING SPECIFIC TO AAAAAA/BBBB IMAGES...
;IF (instrume EQ "BBBB") THEN BEGIN
;   (Edit for other telescopes/instruments)
;ENDIF


; ----------------------------------------------------------------------

; Fit the data in x and y seperately with iterative sigma-clipping and
; construct a "data model" containing the bulk of the signal. Ideally, the
; residuals image contains only shotnoise and the herring-bone pattern.
printf,25,'  constructing data model...'
print,'  constructing data model...'
FLUSH, 25
model = img*0.0

printf,25,'    fitting along columns: yorder=',yorder,' ...',f='(a,i2,a)'
print,'    fitting along columns: yorder=',yorder,' ...',f='(a,i2,a)'
FLUSH, 25
; transpose: swap x and y...
tmodel = TRANSPOSE(model)
timg   = TRANSPOSE(img)

xbuf = INDGEN(nrow)
print,'    ',f='(a,$)'
FOR col=0,ncol-1 DO BEGIN
   ybuf = timg[*,col]

   ; for a stable first-try fit, apply slight smoothing to data...
   mbuf = MEDIAN(ybuf,3,/EVEN)

   ; fix the end-region (virtual overscan)...
   medy = MEDIAN(mbuf,/EVEN)
   mbuf[nrow-2*nos:nrow-1] = medy

   ; first-try cubic spline fit...
   tmpf = SPLFIT(xbuf,mbuf,yorder)

   ; progress indicator ...
   IF (((col+1) MOD 100) EQ 0) THEN BEGIN
      IF ((col+1) LT 1000) THEN BEGIN
         print,(col+1),f='(i3,$)'
      ENDIF ELSE BEGIN
         print,(col+1),f='(i4,$)'
      ENDELSE
   ENDIF ELSE BEGIN   
      IF (((col+1) MOD 25) EQ 0) THEN print,'.',f='(a,$)'
   ENDELSE

   ; iterate fit...
   FOR niter=1,nrej DO BEGIN
      ; reject pixels deviating by more than scut*sigma ...
      rms = STDDEV((ybuf-tmpf))
      ix  = WHERE(ABS(ybuf-tmpf) GT (scut*rms),ixn)
      IF (ixn NE 0) THEN ybuf[ix] = tmpf[ix]
      ; improved fit ...
      tmpf = SPLFIT(xbuf,ybuf,yorder)
   ENDFOR
   tmodel[*,col] = tmpf

ENDFOR
print,' OK'

; transpose image back again...
model = TRANSPOSE(tmodel)
img   = TRANSPOSE(timg)

printf,25,'    fitting along rows:    xorder=',xorder,' ...',f='(a,i2,a)'
print,'    fitting along rows:    xorder=',xorder,' ...',f='(a,i2,a)'
FLUSH, 25
xbuf = indgen(ncol)
print,'    ',f='(a,$)'
FOR row=0,nrow-1 DO BEGIN
   ybuf = img[*,row] - model[*,row]

   ; for a stable first-try fit, apply slight smoothing to data...
   mbuf = MEDIAN(ybuf,3,/EVEN)
   medx = MEDIAN(mbuf,/EVEN)

   ; first-try cubic spline fit...
   tmpf = SPLFIT(xbuf,mbuf,xorder)

   ; progress indicator ...
   IF (((row+1) MOD 100) EQ 0) THEN BEGIN
      IF ((row+1) LT 1000) THEN BEGIN
         print,(row+1),f='(i3,$)'
      ENDIF ELSE BEGIN
         print,(row+1),f='(i4,$)'
      ENDELSE
   ENDIF ELSE BEGIN
      IF (((row+1) MOD 25) EQ 0) THEN print,'.',f='(a,$)'
   ENDELSE

   ; iterate fit...
   FOR niter=1,nrej DO BEGIN
      ; reject pixels deviating by more than scut*sigma ...
      rms = STDDEV((ybuf-tmpf))
      ix  = WHERE(ABS(ybuf-tmpf) GT (scut*rms),ixn)
      IF (ixn NE 0) THEN ybuf[ix] = tmpf[ix]
      ; improved fit ...
      tmpf = SPLFIT(xbuf,ybuf,xorder)
   ENDFOR
   model[*,row] = model[*,row] + tmpf

ENDFOR
print,'  OK'

; end of line-by-line x- and y-fits to the data.


; ----------------------------------------------------------------------

img = TEMPORARY(img) - model
; Note that the mean of array 'img' is now zero.

; Filter out pixels deviating more than 'scut' sigma and add them into
; the data model as well (choose 'scut' sufficiently large that only
; genuine signal (e.g. due to objects, spectral lines, cosmic rays) is
; found, but not the pattern noise; neighbors of rejected pixels may be
; rejected as well depending on the more relaxed sigma criterion 'scnb'.

rms = STDDEV(img)
IF (keyword_set(debug)) THEN print,"      rms(0) = ",rms,f='(a,f10.4)'
FOR niter=1,(nrej-2) DO BEGIN
   rms = STDDEV(img(WHERE(ABS(img) LT (scut*rms))))
   IF (keyword_set(debug)) THEN BEGIN
      print,"      rms(",niter,") = ",rms,f='(a,i1,a,f10.4)'
   ENDIF
ENDFOR

printf,25,'    filtering pixels deviating more than ',scut,$
	'* sigma  [sigma =',rms,'] ...',f='(a,f4.1,a,f7.4,a)'
print,'    filtering pixels deviating more than ',scut,$
	'* sigma  [sigma =',rms,'] ...',f='(a,f4.1,a,f7.4,a)'
FLUSH, 25
ix = WHERE(ABS(img) GT (scut*rms),ixn)
iy = WHERE((ABS(img) GT (scnb*rms) AND ABS(img) LE (scut*rms)),iyn)

IF (ixn NE 0) THEN BEGIN
   ; we found some significantly deviant pixels...

   ; the following lines produce a mask for fast comparison of positions
   ; of the significantly deviant pixels and slightly deviant pixels that
   ; may be neighbors. I fudged the growing radius a bit...
   tmpx = 0.0*img
   tmpx[ix] = ABS(img[ix])
   IF ((ngrow MOD 2) EQ 0) THEN BEGIN
      tmpx = TEMPORARY(SMOOTH(tmpx,(2*ngrow-1)))
   ENDIF ELSE tmpx = TEMPORARY(SMOOTH(tmpx,(2*ngrow+1)))

   ; reject the pixels from the image and place 'em into the model...
   model[ix] = model[ix] + img[ix]
   img[ix]   = 0.0

   printf,25,'      rejected ',ixn,' significantly deviant pixels',f='(a,i6,a)'
   print,'      rejected ',ixn,' significantly deviant pixels',f='(a,i6,a)'
   FLUSH, 25
ENDIF

IF (iyn NE 0 AND ixn NE 0) THEN BEGIN
   ; we found significantly deviant pixels and somewhat deviant pixels;
   ; check whether the latter are neighbors of the former. If so, reject
   ; them too.
   tmpy = 0.0*img
   tmpy[iy] = ABS(img[iy])
   iz = WHERE((tmpx GT 0 AND tmpy GT 0),izn)
   IF (izn NE 0) THEN BEGIN
      model[iz] = model[iz] + img[iz]
      img[iz]   = 0.0
      printf,25,'      rejected ',izn,' neighbors of significantly deviant pixels',$
		f='(a,i6,a)'
      print,'      rejected ',izn,' neighbors of significantly deviant pixels',$
		f='(a,i6,a)'
      FLUSH, 25
   ENDIF
ENDIF


; ----------------------------------------------------------------------

; CONVERT 2-D IMAGE INTO A 1-D TIME SERIES, PADDING THE DEAD-TIME DURING
; THE PARALLEL SHIFT BETWEEN READ-OUT OF INDIVIDUAL PIXEL COLUMNS...
printf,25,'  converting 2-D image to 1-D time-series...'
print,'  converting 2-D image to 1-D time-series...'
FLUSH, 25
nx = ncol + pps
; The non-integer nx prevents phase wandering (TMB).
time_series = FLTARR(nx*nrow)
ds = FLTARR(pps)
FOR row=0,nrow-1 DO BEGIN
   time_series[row*nx] = img[*,row]
   time_series[row*nx+ncol] = ds + MEDIAN(img[*,row],/EVEN)
ENDFOR


; PADD THE TIME-SERIES VECTOR TO AN EVEN NUMBER OF ELEMENTS WITH A SMALL
; PRIME FACTOR, AND PERFORM FORWARD FFT.  Although not as fast as the
; case where the number of elements is a power of 2, this still avoids
; FFT to actually perform a slow FT !  Padding to the next power of 2 is
; impractical, because the number of elements is only a little over 2^20
; pixels for STIS CCD data.
printf,25,'  generating power spectrum (forward FFT) ...'
print,'  generating power spectrum (forward FFT) ...'
FLUSH, 25
IF (keyword_set(debug)) THEN TF0 = SYSTIME(1)
nt = N_ELEMENTS(time_series)
IF ((nt MOD 2L) EQ 1) THEN BEGIN
   ntp = 0
   IF ( nt EQ 1139099 ) THEN BEGIN
      ; 'raw' STIS data...
      ; (1139099+13) = 2^3*3^2*15821     (TMB)
      ; (1139099+37) = 2^6*3*5933        (RAJ: a bit faster [~2.5 x faster])
      ; (1139099+61) = 2^3*3*5*9493      (RAJ: a bit faster [~1.5 x faster])
      ; (1139099+101)= 2^9*5^2*89        (RAJ: very fast!!  [~135 x faster])
      ; (1139099+151)= 2*3*5^3*7^2*31    (RAJ: fastest!!    [~175 x faster])
      ; (1139099+229)= 2^7*3^2*23*43     (RAJ: very fast!!  [~150 x faster])
      ntp = nt + 151
   ENDIF
   IF ( nt EQ 1117277 ) THEN BEGIN
      ; 'flt' STIS data...
      ; (1117277+ 7) = 2^2*3*7*47*283    (TMB)
      ; (1117277+37) = 2*3^5*11^2*19     (RAJ: fastest!!    [~5.5 x faster])
      ; (1117277+151)= 2^2*3*13^2*19*29  (RAJ: a bit faster [~3.5 x faster])
      ; (1117277+379)= 2^3*3^2*19^2*43   (RAJ: a bit faster [~3.4 x faster])
      ntp = nt + 37
   ENDIF
   ;IF ( nt EQ xxxxxx ) THEN BEGIN
   ;  ; adapt to other data...
   ;   ntp = nt + xx
   ;ENDIF
   IF ( ntp EQ 0 ) THEN BEGIN
      print,'WARNING: nt =',nt,': unknown image dimension.'
      print,'	Determine prime number to add to make nt even and optimize'
      print,'	for speed, then edit script AUTOFILET accordingly.'
      RETURN
   ENDIF
ENDIF ELSE BEGIN
   ntp = nt + 6
ENDELSE

t2 = FLTARR(ntp)
t2[0] = time_series

; forward Fourier transform...
ftran = FFT(t2)

IF (keyword_set(debug)) THEN BEGIN
   TF1 = SYSTIME(1)
   print,'    FFT: elapsed time =',(TF1-TF0),' sec',f='(a,f5.2,a)'
ENDIF

; compute frequency array and frequency power spectrum...
freq = FINDGEN(ntp)/(ntp*sst*1.0e-6)
freq[(ntp/2+1):(ntp-1)] = REVERSE(freq[1:(ntp/2-1)])
magn = (ABS(ftran))[0:(ntp/2-1)]	; factor 2 absent in normalization
IF (dc EQ 0) THEN magn[0] = 0		; set power in 0-th element to 0


; save output arrays containing original frequency power spectrum...
frequency  = freq[0:(ntp/2-1)]
omagnitude = magn


; ----------------------------------------------------------------------

; FIND THE PEAK CORRESPONDING TO THE HERRING-BONE PATTERN, AND FIT THE
; PEAK PROFILE. The search interval is restricted by 'freq1' and 'freq2'.
printf,25,'  fitting the peak corresponding to the herring-bone pattern...'
print,'  fitting the peak corresponding to the herring-bone pattern...'
FLUSH, 25
; estimate of width of peak - change for other data
IF (telescop EQ "HST" AND instrume EQ "STIS") THEN psig = 5.5
IF (telescop EQ "TILLINGHAST" AND instrume EQ "FAST") THEN psig = X.X
IF (telescop EQ "AAAAAA" AND instrume EQ "BBBB") THEN psig = X.X

ix = WHERE((frequency GT freq1) AND (frequency LT freq2))
imax = WHERE(omagnitude[ix] EQ MAX(omagnitude[ix]))
peak = omagnitude[ix[imax]]
pctr = frequency[ix[imax]]
params = [peak,pctr,psig,0.0,0.0,0.0]
fparams= [1,1,1,1,0,0]
pfit = CURVEFIT(frequency[ix],omagnitude[ix],frequency[ix]*0+1,params,sigvals,$
	fita=fparams,function_name='rfunct')
coefs = params

ix = WHERE((frequency GT (coefs[1]-12*coefs[2])) AND (frequency LT (coefs[1]+12*coefs[2])))
fparams= [1,1,1,1,1,1]
pfit = CURVEFIT(frequency[ix],omagnitude[ix],frequency[ix]*0+1,coefs,sigvals,$
	fita=fparams,function_name='rfunct')
printf,25,'    Gaussian:  ampl=',coefs[0],',  cntr=',coefs[1],$
	' Hz,  sigm=',coefs[2],' Hz',f='(a,f11.4,a,f11.4,a,f11.4,a)'
printf,25,'    Backgrnd:  mean=',coefs[3],',  tilt=',coefs[4],$
	',     curv=',coefs[5],'   ',f='(a,f11.6,a,f11.8,a,f11.8,a)'
print,'    Gaussian:  ampl=',coefs[0],',  cntr=',coefs[1],$
	' Hz,  sigm=',coefs[2],' Hz',f='(a,f11.4,a,f11.4,a,f11.4,a)'
print,'    Backgrnd:  mean=',coefs[3],',  tilt=',coefs[4],$
	'      curv=',coefs[5],'   ',f='(a,f11.6,a,f11.8,a,f11.8,a)'
FLUSH, 25

; ----------------------------------------------------------------------

; TEST FOR ANOMALOUS PATTERN PARAMETERS (frequency,amplitude,drift-width)
; We will use the results from an earlier trending analysis (Jansen et al.
; 2005) to decide whether or not a pattern is anomalous.  That study used  
; a reference date of 2002.25 (mjd=52367) and median frequency 16117.8 Hz,
; and a quadratic polynomial fit to describe the trend of the median
; frequency with time:
mjd0 = 52367.0
frq0 = 16117.8
fdt  = (tmjd - mjd0)/365.25
rfrq = coefs[1] + 351.7*fdt - 33.6*fdt^2
dwdt = coefs[2]
pamp = coefs[0]*FLOAT(gain)

; set frequency, drift-width and amplitude warning flags...
fflag = 0
wflag = 0
aflag = 0
IF ( ABS(rfrq-frq0) GT 1400 ) THEN fflag = 1
IF ((ABS(dwdt-9.5) GT 6.5) OR (dwdt LT 2.0)) THEN wflag = 1
IF ((ABS(pamp-0.12) GT 0.10) OR (pamp LT 0.01)) THEN aflag = 1

; if any of the parameters is anomalous:
IF ((fflag GT 0) OR (aflag GT 0) OR (wflag GT 0)) THEN BEGIN
   printf,25,'    WARNING: anomalous pattern in '+image+'['+STRTRIM(exten,2)+']'+'  VERIFIED=NO  GENUINE=UNKNOWN',f='(a)'
   print,'    WARNING: anomalous pattern in '+image+'['+STRTRIM(exten,2)+']'
   FLUSH, 25
   IF (fflag GT 0) THEN BEGIN
      IF ( (rfrq-frq0) LT 0 ) THEN BEGIN
         printf,25,'    WARNING:   frequency < ',(frq0-1400),' Hz',f='(a,f7.1,a)'
         FLUSH, 25
      ENDIF ELSE BEGIN
         printf,25,'    WARNING:   frequency > ',(frq0+1400),' Hz',f='(a,f7.1,a)'
         FLUSH, 25
      ENDELSE
   ENDIF
   IF (wflag GT 0) THEN BEGIN
      IF ( dwdt LT 2.0 ) THEN BEGIN
         printf,25,'    WARNING:   driftwidth < 2.0 Hz  (too small)',f='(a)'
         FLUSH, 25
      ENDIF ELSE BEGIN
         IF ( dwdt GT 16.0 ) THEN BEGIN
            printf,25,'    WARNING:   driftwidth > 16.0 Hz (too large?)',f='(a)'
            FLUSH, 25
         ENDIF ELSE BEGIN
            printf,25,'    WARNING:   driftwidth < 3.0 Hz  (too small?)',f='(a)'
            FLUSH, 25
         ENDELSE
      ENDELSE
   ENDIF
   IF (aflag GT 0) THEN BEGIN
      IF ( pamp LT 0.01 ) THEN BEGIN
         printf,25,'    WARNING:   amplitude < 0.01  (undetected?)',f='(a)'
         FLUSH, 25
      ENDIF ELSE BEGIN
         IF ( pamp GT 0.22 ) THEN BEGIN
            printf,25,'    WARNING:   amplitude > 0.22  (too large?)',f='(a)'
            FLUSH, 25
         ENDIF ELSE BEGIN
            printf,25,'    WARNING:   amplitude < 0.02  (too small?)',f='(a)'
            FLUSH, 25
         ENDELSE
      ENDELSE
   ENDIF
ENDIF



; ----------------------------------------------------------------------

;"WIPE" FREQUENCIES RESPONSIBLE FOR HERRING-BONE PATTERN...
; NOTE: this is not identical to the original 'wipe' function of TMB: here,
;       all signal in the frequency interval associated with the herring-bone
;       pattern will be substituted by white noise that matches the noise
;       outside this frequency interval, rather than suppressed (scaled) with
;       a (top-hat or tapered) window function.
printf,25,'  wiping frequencies responsible for herring-bone pattern...'
print,'  wiping frequencies responsible for herring-bone pattern...'
FLUSH, 25
xlim = [(coefs[1]-3*coefs[2]),(coefs[1]+3*coefs[2])]
xliml= [(coefs[1]-7*coefs[2]),(coefs[1]-4*coefs[2])]
xlimh= [(coefs[1]+4*coefs[2]),(coefs[1]+7*coefs[2])]
ix   = WHERE((freq GE xlim[0]) AND (freq LE xlim[1]))
ixout= WHERE(((freq GE xliml[0]) AND (freq LE xliml[1])) $
		OR ((freq GE xlimh[0]) AND (freq LE xlimh[1])))
ampl = STDDEV(magn[ixout])
seed = SYSTIME(1)
filt = ftran
; function RANDOMN returns normally distributed pseudo-random numbers with a
; mean of zero and a standard deviation of one.
filt[ix] = ampl*RANDOMN(seed, N_ELEMENTS(ix)) + coefs[3]

; save output array containing filtered frequency power spectrum...
fmagnitude = (ABS(filt))[0:(ntp/2-1)]


; ----------------------------------------------------------------------

; write original and filtered power spectra and fit to output table...
printf,25,'    saving original and filtered spectra to "'+fpspec+'"...'
print,'    saving original and filtered spectra to "'+fpspec+'"...'
FLUSH, 25
; write header of table with frequency power spectra...
OPENW,26,fpspec
printf,26,'# freq omagn fmagn coefs',f='(a)'
printf,26,'# ',f='(a)'
printf,26,'# ',progver,f='(a,a)'
printf,26,'#   image  ="',image,'"',f='(a,a,a)'
printf,26,'#   result ="',result,'"',f='(a,a,a)'
printf,26,'# ',f='(a)'
printf,26,'# Fit:  F(f) = c0*EXP{-0.5*(f/c2)^2} + c3 + c4*f + c5*f^2',f='(a)'
printf,26,'#       where f = (freq-c1)',f='(a)'
printf,26,'# ',f='(a)'
printf,26,'#   f (Hz)      omagn      fmagn          c     ',f='(a)'
printf,26,'#===============================================',f='(a)'

; write frequency power spectra and fit coefficients...
ix = WHERE(frequency GT (coefs(1)-10*coefs(2)) AND frequency LT (coefs(1)+10*coefs(2)))
FOR nn=0,N_ELEMENTS(ix)-1 DO BEGIN
   IF ( nn LE 5 ) THEN BEGIN
      printf,26,frequency[ix[nn]],omagnitude[ix[nn]],fmagnitude[ix[nn]],$
	coefs[nn],f='(f11.4,f11.6,f11.6,f15.8)'
   ENDIF ELSE BEGIN
      printf,26,frequency[ix[nn]],omagnitude[ix[nn]],fmagnitude[ix[nn]],$
	f='(f11.4,f11.6,f11.6)'
   ENDELSE
ENDFOR

CLOSE,26,/FORCE


; ----------------------------------------------------------------------

; PERFORM INVERSE FFT AND GENERATE RESULT (CLEANED) AND PATTERN IMAGES...
printf,25,'  generating result and herring-bone pattern image (inverse FFT)...'
print,'  generating result and herring-bone pattern image (inverse FFT)...'
FLUSH, 25
btran = FFT(filt,/INVERSE)
ttf = (FLOAT(btran))[0:nt]
; generate output images...
outimg = FLTARR(ncol,nrow)
herimg = FLTARR(ncol,nrow)
FOR row=0,nrow-1 DO BEGIN
  outimg[*,row] = model[*,row] + ttf[row*nx:(row*nx+ncol-1)]
  herimg[*,row] = img[*,row] - ttf[row*nx:(row*nx+ncol-1)]
ENDFOR
IF (type EQ "flt" OR type EQ "trimmed") THEN BEGIN
   outimg = outimg[nos:(nc0-nos-1),*]
   herimg = herimg[nos:(nc0-nos-1),*]
ENDIF

; RESTORE ORIGINAL IMAGE ORIENTATION...
IF (instrume EQ "STIS") THEN BEGIN
   IF (amp EQ "B") THEN outimg = REVERSE(outimg,1)
   IF (amp EQ "C") THEN outimg = REVERSE(outimg,2)
   IF (amp EQ "D") THEN outimg = ROTATE(outimg,2)
   IF (amp EQ "B") THEN herimg = REVERSE(herimg,1)
   IF (amp EQ "C") THEN herimg = REVERSE(herimg,2)
   IF (amp EQ "D") THEN herimg = ROTATE(herimg,2)
ENDIF


; ----------------------------------------------------------------------

; WRITE OUTPUT FITS IMAGE(S)...
T1a = SYSTIME(0)
GET_DATE,dte,/TIMETAG

; First update primary header to conform to the FITS standard...
;; note, that there a problem with breaking the arguments of SXADDPAR over
;; multiple lines using a '$' (RAJ - Aug 8, 2002)

SXADDPAR,primehead,'SIMPLE','T',' FITS standard'
SXADDPAR,primehead,'BITPIX',-32,' Bits per data value'
SXADDPAR,primehead,'NAXIS',2,' Number of data axes'
SXADDPAR,primehead,'NAXIS1',ncol,' Length of axis 1',BEFORE='EXTEND'
SXADDPAR,primehead,'NAXIS2',nrow,' Length of axis 2',BEFORE='EXTEND'
SXADDPAR,primehead,'EXTEND','F',' File may contain standard extensions'
SXDELPAR,primehead,'NEXTEND'
SXDELPAR,primehead,'GROUPS'
SXADDPAR,primehead,'ORIGIN',progver,' FITS file originator',BEFORE='DATE'
SXADDPAR,primehead,'DATE',dte,' Date this file was written'
SXADDPAR,primehead,'FILENAME',STRMID(result,STRPOS(result,'/',/REVERSE_SEARCH)+1),' Name of FITS file'
SXADDPAR,primehead,'HPATCORR','COMPLETE',' Remove herring-bone pattern noise',BEFORE='DQICORR'
SXDELPAR,primehead,'END'
SXDELPAR,primehead,'HISTORY'
SXDELPAR,primehead,'COMMENT'

SXDELPAR,extenhead,'XTENSION'
SXDELPAR,extenhead,'BITPIX'
SXDELPAR,extenhead,'NAXIS'
SXDELPAR,extenhead,'NAXIS1'
SXDELPAR,extenhead,'NAXIS2'
SXDELPAR,extenhead,'PCOUNT'
SXDELPAR,extenhead,'GCOUNT'
SXDELPAR,extenhead,'INHERIT'
SXDELPAR,extenhead,'EXTNAME'
SXDELPAR,extenhead,'EXTVER'
SXDELPAR,extenhead,'DATE'
SXDELPAR,extenhead,'ORIGIN'

newhead = [primehead,extenhead]

SXADDPAR,newhead,'HISTORY',' '
SXADDPAR,newhead,'HISTORY',' HPATCORR Started '+T0a
SXADDPAR,newhead,'HISTORY','  '+progver
SXADDPAR,newhead,'HISTORY','   image="'+STRMID(image,STRPOS(image,'/',/REVERSE_SEARCH)+1)+'['+STRTRIM(STRING(exten),2)+']'+'", interval='+STRTRIM(STRING(freq1),2)+'--'+STRTRIM(STRING(freq2),2)+' Hz'
SXADDPAR,newhead,'HISTORY','   xorder='+STRTRIM(STRING(xorder),2)+', yorder='+STRTRIM(STRING(yorder),2)+', scut='+STRMID(STRTRIM(STRING(scut),2),0,4)+', scnb='+STRMID(STRTRIM(STRING(scnb),2),0,4)+', nrej='+STRTRIM(STRING(nrej),2)+', ngrow='+STRTRIM(STRING(ngrow),2)
SXADDPAR,newhead,'HISTORY','  Pattern characteristics:'
SXADDPAR,newhead,'HISTORY','   freq='+STRTRIM(STRING(coefs[1]),2)+' Hz, width='+STRTRIM(STRING(coefs[2]),2)+' Hz, pmax='+STRTRIM(STRING(coefs[0]),2)
SXADDPAR,newhead,'HISTORY',' HPATCORR Finished '+T1a

printf,25,'  writing FITS image "',result,'"...'
print,'  writing FITS image "',result,'"...'
FLUSH, 25
WRITEFITS,result,outimg,newhead


; WRITE FITS IMAGES FOR INTERMEDIATE STAGES (DEBUGGING)...
IF (keyword_set(debug)) THEN BEGIN
   printf,25,'  exporting data model, residuals, and pattern-noise images...'
   print,'  exporting data model, residuals, and pattern-noise images...'
   FLUSH, 25
   model = rotate(model,2)
   img   = rotate(img,2)
   tmpname = STRMID(result,0,STRPOS(result,'_'))+'_res_'+STRTRIM(STRING(exten),2)+'.fits'
   ;tmpname = STRMID(result,0,STRPOS(result,'.',/REVERSE_SEARCH))+"_res.fits"
   print,'    writing FITS image "'+tmpname+'"...'
   SXADDPAR,newhead,'FILENAME',STRMID(tmpname,STRPOS(tmpname,'/',/REVERSE_SEARCH)+1),' Name of FITS file'
   WRITEFITS,tmpname,img,newhead
   tmpname = STRMID(result,0,STRPOS(result,'_'))+'_mod_'+STRTRIM(STRING(exten),2)+'.fits'
   ;tmpname = STRMID(result,0,STRPOS(result,'.',/REVERSE_SEARCH))+"_mod.fits"
   print,'    writing FITS image "'+tmpname+'"...'
   SXADDPAR,newhead,'FILENAME',STRMID(tmpname,STRPOS(tmpname,'/',/REVERSE_SEARCH)+1),' Name of FITS file'
   WRITEFITS,tmpname,model,newhead
   tmpname = STRMID(result,0,STRPOS(result,'_'))+'_rwh_'+STRTRIM(STRING(exten),2)+'.fits'
   ;tmpname = STRMID(result,0,STRPOS(result,'.',/REVERSE_SEARCH))+"_rwh.fits"
   print,'    writing FITS image "'+tmpname+'"...'
   SXADDPAR,newhead,'FILENAME',STRMID(tmpname,STRPOS(tmpname,'/',/REVERSE_SEARCH)+1),' Name of FITS file'
   WRITEFITS,tmpname,herimg,newhead

ENDIF


T1 = SYSTIME(1)

; Closing message...
printf,25,'  AUTOFILET: finished (elapsed time: ',(T1-T0),' sec).',f='(a,f7.2,a)'
print,'  AUTOFILET: finished (elapsed time: ',(T1-T0),' sec).',f='(a,f7.2,a)'
FLUSH, 25
CLOSE,25,/ALL,/FORCE

RETURN
END

;-------------------------------  R A J  -------------------------------
