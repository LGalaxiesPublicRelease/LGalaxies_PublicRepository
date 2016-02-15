;+
; NAME:
;   bvalu2d
;
; PURPOSE:
;   Evaluate a bspline fit over 2 dependent variables
;
; CALLING SEQUENCE:
;   
;    z = bvalu2d(x, y, bkpt, coeff, ideriv=ideriv)
;
; INPUTS:
;   x          - vector of x positions to evaluate
;   y          - vector of y positions to evaluate
;   bkpt       - Breakpoint vector returned by efc
;   coeff      - B-spline coefficients calculated by efc 
;		   2d in this case [npoly, nbkpt-nord]
;
; OPTIONAL KEYWORDS:
;
;   ideriv     - Derivative to evaluate at x (default 0)
;
; OUTPUTS:
;   z          - Evaluations corresponding to x positions
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES CALLED:
;   slatec_bvalu
;
;
; REVISION HISTORY:
;   13-Mar-2000  Written by Scott Burles, FNAL
;-
;------------------------------------------------------------------------------
;
;	Coeff actually contains a nord x ncoeff numbers
;       containing 2d information
;


function bvalu2d, x, y, fullbkpt, coeff, ideriv=ideriv

    npoly = (size(coeff))[0]
    ncoeff = (size(coeff))[1]

    z = x*0.0

    for i=0,npoly -1 do begin

      z = z + slatec_bvalu(x, fullbkpt, coeff[i,*], ideriv=ideriv) * y^i
    endfor

    return, z
end
