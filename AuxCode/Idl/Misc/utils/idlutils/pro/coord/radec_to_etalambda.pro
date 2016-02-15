;+
; NAME:
;   radec_to_etalambda
; PURPOSE:
;   convert from RA, Dec to eta, lambda (SDSS survey coordinates)
; CALLING SEQUENCE:
;   radec_to_etalambda, ra,dec,eta,lambda,stripenumber=stripenumber
; INPUTS:
;   ra      RA (deg), J2000
;   dec     Dec (deg), J2000
; OUTPUTS:
;   eta     SDSS survey coordinate eta (deg)
;   lambda  SDSS survey coordinate lambda (deg)
; OPTIONAL OUTPUTS:
;   stripenumber   SDSS survey stripe number (integer)
; BUGS:
;   Location of the survey center is hard-wired (in etalambda_to_radec.pro);
;     it should be read from astrotools.
; REVISION HISTORY:
;   2001-Jul-21  written by Hogg (NYU)
;   2002-Oct-04  modified by Blanton (NYU)
;-
pro radec_to_etalambda, ra,dec,eta,lambda,stripenumber=stripenumber

if (n_params() NE 4) then begin
  print, 'Syntax - radec_to_etalambda, ra, dec, eta, lambda [, stripenumber=]'
  return
endif
if (n_elements(ra) NE n_elements(dec)) then $
  message, 'Number of elements in RA and DEC must agree'

r2d= 180.0D/(!DPI)
d2r= 1.D/r2d

; find center
etalambda_to_radec, 0.0D,0.0D,racen,deccen

; basic transformation
eta= r2d*atan(sin(dec*d2r),cos(dec*d2r)*cos((ra-racen)*d2r))-deccen
stripenumber= floor((eta+58.75D)/2.5D)
lambda= r2d*asin(cos(dec*d2r)*sin((ra-racen)*d2r))

; all the sign-flipping crap
bad= where(eta LT -90.0D OR eta GT 90.0D,nbad)
if nbad GT 0 then begin
    eta[bad]= eta[bad]+180.0D
    lambda[bad]= (180.0D)-lambda[bad]
endif

bad=where(eta gt 180.D,nbad)
if(nbad gt 0) then eta[bad]=eta[bad]-360.D

bad=where(lambda gt 180.D,nbad)
if(nbad gt 0) then lambda[bad]=lambda[bad]-360.D

bad=where(stripenumber lt 0.D,nbad)
if(nbad gt 0) then stripenumber[bad]=stripenumber[bad]+144L

end

