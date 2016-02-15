;+
; NAME:
;   psf_norm
;
; PURPOSE:
;   Sum the counts inside a box of radius rad centered on stamp
;
; CALLING SEQUENCE:
;   norm = psf_norm(stamp, rad)
;
; INPUTS:
;   stamp    - postage stamp of star
;   rad      - radius of normalization box
;
; OUTPUTS:
;   norm     - total of counts in centered box
;
; EXAMPLES:
;   norm = psf_norm(stampcen, par.cenrad)
;
; COMMENTS:
;   Sum the counts inside a box of radius rad centered on stamp.
;   This is thought to be more stable than using the peak value. 
;   To use peak value, simply set rad to zero.
;
; REVISION HISTORY:
;   2006-May-25   Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function psf_norm, stamp, rad

  if NOT keyword_set(rad) then message, 'must set rad'
  sz = size(stamp, /dimen)
  cen = (sz-1)/2

  if rad GT min(cen) then message, 'rad too big'

  subimg = stamp[cen[0]-rad:cen[0]+rad, cen[1]-rad:cen[1]+rad]
  norm   = total(subimg)

  return, norm
end
