;+
; NAME:
;   frac_profmean
; PURPOSE:
;   Get the radius of a certain light fraction, given the total light radius 
; CALLING SEQUENCE
;   radfrac= frac_profmean(frac, nprof, profmean, radtot $
;                          [,profradius=profradius] )
; INPUTS:
;   frac        - fraction desired
;   nprof       - number of annuli to trust
;   profmean    - [nrad,ncent] annular fluxes
;   radtot      - total light radius 
; OPTIONAL INPUTS:
;   profradius  - [nrad] defining profile, in pixels; default to PHOTO aps
; KEYWORDS:
; OUTPUTS:
;   radfrac     - radius containing frac fraction of the light in radtot
; COMMENTS:
; DEPENDENCIES:
;   idlutils
; BUGS:
;   no indication of whether fit is sensible
; REVISION HISTORY:
;   2003-09-15  Written - Blanton
;-
function frac_profmean_func, rad

common frac_profmean_com, frac, nprof, profmean, fluxtot, profradius

interp_profmean,nprof,profmean,rad,flux
return,fluxtot*frac-flux

end
;
function frac_profmean, in_frac, in_nprof, in_profmean, radtot, $
                        profradius=in_profradius

common frac_profmean_com 

; set defaults
if NOT keyword_set(in_profradius) then $
  in_profradius= [  0.564190,   1.692569,   2.585442,   4.406462, $
                 7.506054,  11.576202,  18.584032,  28.551561, $ 
                 45.503910,  70.510155, 110.530769, 172.493530, $
                 269.519104, 420.510529, 652.500061]

frac=in_frac
nprof=in_nprof
profmean=in_profmean
profradius=in_profradius

if(nprof eq 0) then return,0.

interp_profmean,nprof,profmean,radtot,fluxtot
if(fluxtot eq !VALUES.F_INFINITY) then $
  return,radtot
radfrac=zbrent(profradius[0]*0.01,radtot,func_name='frac_profmean_func')

return, radfrac

end
