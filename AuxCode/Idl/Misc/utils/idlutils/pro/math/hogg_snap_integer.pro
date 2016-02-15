;+
; NAME:
;   hogg_snap_integer
; PURPOSE:
;   Randomly snap floats or doubles to nearest integers up and down so
;     that the mean is correct on average.
; CALLING SEQUENCE:
;   ii= hogg_snap_integer(seed,xx)
; INPUTS:
;   seed - long seed for randomu
;   xx   - array of floating point numbers
; OUTPUTS:
;   ii   - array of integers
; BUGS:
;   Relies on "true" being 1.
; REVISION HISTORY:
;   2003-08-13  written - Hogg
;-
function hogg_snap_integer, seed,xx
dimens= size(xx,/dimens)
ii= long(djs_floor(xx))
dx= xx-double(ii)
if dimens[0] EQ 0 then $
  ii= ii+(randomu(seed) LT dx) $
else $
  ii= ii+(randomu(seed,dimens) LT dx)
return, ii
end
