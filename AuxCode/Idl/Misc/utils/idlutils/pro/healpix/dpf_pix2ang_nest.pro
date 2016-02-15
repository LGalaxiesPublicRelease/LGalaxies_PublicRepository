;+
; NAME:
;   dpf_pix2ang_nest
;
; PURPOSE:
;   Compute coordinates for HEALPix pixel numbers, nested order
;
; CALLING SEQUENCE:
;   dpf_pix2ang_nest, nside, ipix, theta, phi, double=double
;
; INPUTS:
;   nside   - healpix nside
;   ipix    - pixel numbers
; 
; OUTPUTS:
;   theta   - angle from north pole [radians]
;   phi     - longitude angle [radians]
;
; EXAMPLES:
;   
; COMMENTS:
;   Calling syntax is same as pix2ang_nest and agrees to machine 
;     precision. 
;   This routine has somewhat better performance when called for the
;     full sky than pix2ang_nest. 
;   
; REVISION HISTORY:
;   2004-Jun-08  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
pro dpf_pix2ang_nest, nside, ipix, theta, phi, double=double

  t1 = systime(1)

  nside  = long(nside)
  npix   = 12L*nside*nside
  nlevel = round(alog(nside)/alog(2))
  if 2L^nlevel NE nside then message, 'bad nside value - must be power of 2.'

  dpf_nest2ring, nside, ipix, ipring
  dpf_pix2ang_ring, nside, ipring, theta, phi, double=double

; -------- number in each nest
;  print, 'time', systime(1)-t1

  return
end
