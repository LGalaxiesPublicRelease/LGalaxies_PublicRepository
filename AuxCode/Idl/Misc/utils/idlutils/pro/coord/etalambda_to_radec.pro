;+
; NAME:
;   etalambda_to_radec
; PURPOSE:
;   convert from eta, lambda (SDSS survey coordinates) to RA, Dec
; CALLING SEQUENCE:
;   etalambda_to_radec, eta,lambda,ra,dec
; INPUTS:
;   eta     SDSS survey coordinate eta (deg)
;   lambda  SDSS survey coordinate lambda (deg)
; OPTIONAL OUTPUTS:
;   ra      RA (deg), J2000
;   dec     Dec (deg), J2000
; BUGS:
;   Location of the survey center is hard-wired, not read from astrotools.
; REVISION HISTORY:
;   2001-Jul-21  written by Hogg (NYU)
;   2002-Oct-04  modified by Blanton (NYU)
;-
pro etalambda_to_radec, eta,lambda,ra,dec

if (n_params() NE 4) then $
  message, 'Wrong number of parameters'
if (n_elements(ra) NE n_elements(dec)) then $
  message, 'Number of elements in RA and DEC must agree'

racen= 185.0D                   ; deg
deccen= 32.5D                   ; deg
r2d= 180.0D/(!DPI)
d2r= 1.D/r2d

coslambda=cos(lambda*d2r)
dec= r2d*asin(coslambda*sin((eta+deccen)*d2r))
ra= r2d*atan(sin(lambda*d2r),coslambda*cos((eta+deccen)*d2r))+racen

end
