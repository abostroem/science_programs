; $Id: //depot/idl/IDL_63_RELEASE/idldir/lib/rfunct.pro#1 $
;
; Copyright (c) 1982-2006, Research Systems, Inc.  All rights reserved.
;	Unauthorized reproduction prohibited. Modified (c)2007, R.A. Jansen.
;

PRO	RFUNCT,X,A,F,PDER
;+
; NAME:
;	RFUNCT
;
; PURPOSE:
;	Evaluate the sum of a Gaussian and a 2nd-order polynomial
;	and optionally return the value of its partial derivatives.
;	Normally, this function is used by CURVEFIT to fit the
;	sum of a line and a varying background to actual data.
;
; CATEGORY:
;	E2 - Curve and surface fitting.
;
; CALLING SEQUENCE:
;	RFUNCT, X, A, F [, Pder]
;
; INPUTS:
;	X:	The values of the independent variable.
;	A:	The parameters of the equation described in PROCEDURE below.
;
; OUTPUTS:
;	F:	The value of the function at each X(i).
;
; OPTIONAL OUTPUT PARAMETERS:
;	Pder:	An array of the size (N_ELEMENTS(X),6) that contains the
;		partial derivatives.  Pder(i,j) represents the derivative
;		at the i'th point with respect to j'th parameter.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	None.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	F = A(0)*EXP(-0.5*(Z/A(2))^2) + A(3) + A(4)*Z + A(5)*Z^2
;	Z = (X-A(1))
;
; MODIFICATION HISTORY:
;	WRITTEN, DMS, RSI, SEPT, 1982.
;	Modified, DMS, Oct 1990.  Avoids divide by 0 if A(2) is 0.
;	Modified, RAJ, Aug 2007.  Changed Z from (X-A(1))/A(2) to (X-A(1))
;		and also use Z to fit 1st and 2nd order polynomial terms.
;-
	ON_ERROR,2                        ;Return to caller if an error occurs

	;GET Z
	Z = (X-A[1])
	;GAUSSIAN PART, IGNORE SMALL TERMS
	if a[2] ne 0.0 then begin
	   EZ = EXP(-0.5*(Z/A[2])^2)*(ABS(Z/A[2]) LE 7.)
	endif else begin
	   EZ = EXP(-0.5*10.^2)*(ABS(Z/A[2]) LE 7.)
	endelse

	;FUNCTIONS.
	F = A[0]*EZ + A[3] + A[4]*Z + A[5]*Z^2

	IF N_PARAMS(0) LE 3 THEN RETURN ;NEED PARTIAL?
;
	PDER = FLTARR(N_ELEMENTS(X),6) ;YES, MAKE ARRAY.
	;COMPUTE PARTIALS
	PDER[0,0] = EZ
	if a[2] ne 0. then PDER[0,1] = A[0] * EZ * Z/(A[2])^2
	PDER[0,2] = PDER[*,1] * Z
	PDER[*,3] = 1.
	PDER[0,4] = Z
	PDER[0,5] = Z^2
	RETURN
END

