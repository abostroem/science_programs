function splfit,x,y,nodes,xs,ys
;+
;                       splfit
;
; driver routine for splinefit to compute least squares spline
; with equally spaced nodes.
;
; CALLING SEQUENCE:
;       fit = splfit(x,y,nodes,xs,ys)
;
; INPUTS:
;       x - x data vector
;       y - y data vector
;       nodes - number of spline nodes
;
; OUTPUT:
;       the least square spline, evaluated at points x, is returned
;       as the function value.
;
; OPTIONAL OUTPUTS:
;       xs - x positions of the least squares spline
;       ys - y positions of the least sqaures spline
;
; HISTORY:
;       version 1  D. Lindler   May 1989
;       version 2  DJL/ACC      1/22/92         Moved to IDL Version 2
;       14-JUL-95 Deutsch renamed 'sfit' to 'splfit' since IDL now includes
;               a routine called sfit (surface fitting)
;       version 2a R.A. Jansen  Aug 22 2007     Slight reformatting to improve
;               readability of code.
;-
;---------------------------------------------------------------------------
IF ( n_params(0) LT 3 ) THEN nodes=10

; '+1' in definition of 'ys' because in 'splinefit' DELY = 0.01*ys(i) > 0 ...
IF ( nodes GT 0 ) THEN BEGIN
        xs = FINDGEN(nodes)*(MAX(x)-MIN(x))/(nodes-1) + MIN(x)
        ys = TOTAL(y)/n_elements(y) + xs*0 + 1
END

; uniform weights in fit...
w = x*0 + 1

return,splinefit(x,y,w,xs,ys,sig,0.01)

END
