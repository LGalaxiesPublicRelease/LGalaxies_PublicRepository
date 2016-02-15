;+
; NAME:
;   hogg_tp_shift
; PURPOSE:
;   Shift tangent point on the sphere (CRVAL, in RA, Dec units),
;   adjusting simultaneously the CD matrix (to deal with coordinate
;   issues near the celestial poles).  The idea is to shift the WCS
;   without substantially rotating the tangent-plane coordinates, even
;   when near the poles.
; INPUTS:
;   astr   - astrometry structure
;   crval  - new crval (tangent point on the sphere) to insert
; OUTPUTS:
;          - new astrometry structure
; REVISION HISTORY:
;   2005-08-21  started - Hogg (NYU)
;-
function hogg_tp_shift, astr,crval
ra_old=  astr.crval[0]
ra_new=  crval[0]
dec_new= crval[1]
deltaorient= (ra_new-ra_old)*!DPI/180D0*sin(dec_new*!DPI/180D0)
newastr= astr
newastr.crval= crval
newastr.cd= astr.cd#[[ cos(deltaorient), sin(deltaorient)], $
                     [-sin(deltaorient), cos(deltaorient)]]
return, newastr
end
