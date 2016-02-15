;+
; NAME:
;   heal_smooth
;
; PURPOSE:
;   Smooth ring-ordered healpix maps with spherical harmonic convolution
;
; CALLING SEQUENCE:
;   smooth_map = heal_smooth(map, fwhm_arcmin, nside=, alm=, lmax= )
;
; INPUTS:
;   map         - healpix map
;   fwhm_arcmin - fwhm of gaussian smoothing kernel (arcmin)
;   nside       - set if map not set and alms passed
;   alm         - set if alms already known (saves time)
;   lmax        - maximum l value to use
;
; OPTIONAL INPUTS:
;   
; KEYWORDS:
;   
; OUTPUTS:
;   
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   Input maps must be ring ordered - use nest2ring otherwise
;   
; REVISION HISTORY:
;   2003-Mar-14  Written by Douglas Finkbeiner, Princeton
;   2004-Aug-16  Avoid floating underflow in beam transform - DPF
;
;----------------------------------------------------------------------
function heal_smooth, map, fwhm_arcmin, nside=nside, alm=alm_in, lmax=lmax

  if n_elements(fwhm_arcmin) eq 0 then message, 'must set fwhm_arcmin'
  if fwhm_arcmin le 0 then return, map

; -------- if we pass a map, compute alm transform
  if keyword_set(map) then begin 
     nside_in = long(sqrt(n_elements(map)/12))
     if NOT keyword_set(lmax) then lmax = 2*nside_in
     alm = healpix2alm(map, lmax=lmax)
     if NOT keyword_set(nside) then nside = nside_in
  endif else begin 
; -------- if we pass alm array, simply determine lmax
     if NOT keyword_set(alm_in) then message, 'must set either map or alm'
     if NOT keyword_set(lmax) then lmax = (size(alm_in, /dim))[0]
     if NOT keyword_set(nside) then message, 'must set nside when passing alm'
     alm = alm_in
  endelse

  if nside GT 8192 then message, 'routines cannot handle maps this big'

; -------- sigma in radians
  sigma =  fwhm_arcmin *!dpi/180.d /60.d / sqrt(8.d *alog(2.))

; -------- beam transform
  l   = dindgen(lmax)
  exponent = -0.5d * l*(l+1) * sigma^2
  smallexp = where(exponent LT -300, nsmall)
  bl  = exp(exponent >  (-300))
  if nsmall GT 0 then bl[nsmall] = 0

;  plot, l, bl, title='beam transform'

; -------- apply smoothing to alm array
  lind = where(bl gt 0.001, nl)
  nl = nl < (nside * 4)
  if nl gt 0 then begin 
     for ll=0, nl-1 do alm[ll, *] = alm[ll, *]*bl[ll]
  endif 

  print, 'NL = ', nl
  smooth_map = alm2healpix(nside, alm, lmax=nl-1)

  return, smooth_map
end
