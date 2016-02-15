;+
; NAME:
;   psf_validate
;
; PURPOSE:
;   sanity checks on PSF fit
;
; CALLING SEQUENCE:
;   psf_validate, pstr
;
; INPUTS:
;   pstr   - psf structure (see psf_fit_coeffs
;
; OUTPUTS:
;   junk to screen
;
; OPTIONAL OUTPUTS:
;   
; RESTRICTIONS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   
; REVISION HISTORY:
;   2006-Jun-06   Written by Douglas Finkbeiner, Princeton (6/6/6)
;
;----------------------------------------------------------------------
pro psf_validate, pstr

  psfs = psf_eval(pstr.x, pstr.y, pstr.coeff, pstr.cenrad)

  npix = (size(psfs, /dimen))[0]
  midind = (npix*npix-1)/2

  mask = bytarr(npix, npix)
  mask[pstr.boxrad-pstr.fitrad:pstr.boxrad+pstr.fitrad, $
       pstr.boxrad-pstr.fitrad:pstr.boxrad+pstr.fitrad] = 1B
  mask[pstr.boxrad-pstr.fitrad+1:pstr.boxrad+pstr.fitrad-1, $
       pstr.boxrad-pstr.fitrad+1:pstr.boxrad+pstr.fitrad-1] = 0B
  ringind = where(mask)

  for i=0L, pstr.nstar-1 do begin 

; -------- check central pixel is highest
     psf = psfs[*, *, i]
     junk = max(psf, maxind)
     if maxind NE midind then splog, 'center pixel not highest!'
     noise = stdev(psf[ringind])/psf[midind]
     if noise GT .02 then print, i, noise

  endfor

  return
end
