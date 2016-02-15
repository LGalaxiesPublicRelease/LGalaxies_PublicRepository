;+
; NAME:
;   map_nest2ring
;
; PURPOSE:
;   Convert a full-sky map in nested order to ring order
;
; CALLING SEQUENCE:
;   ringmap=map_nest2ring(nestmap)
;
; INPUTS:
;   nestmap  - healpix map in nested order
;
; OUTPUTS:
;   ringmap  - healpix map in ring order
;
; COMMENTS:
;   
; REVISION HISTORY:
;   2003-Dec-04  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function map_nest2ring, nestmap

  npix = n_elements(nestmap)
  nside = round(sqrt(npix/12))
       
  dpf_nest2ring, nside, lindgen(12*nside*nside), ipring

  ringmap = make_array(npix, type=size(nestmap, /type), /nozero)
  ringmap[ipring] = nestmap

  return, ringmap
end
