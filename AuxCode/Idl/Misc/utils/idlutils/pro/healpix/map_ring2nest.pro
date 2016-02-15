;+
; NAME:
;   map_ring2nest
;
; PURPOSE:
;   Convert a full-sky map in ring order to nested order
;
; CALLING SEQUENCE:
;   nestmap=map_ring2nest(ringmap)
;
; INPUTS:
;   ringmap  - healpix map in ring order
;
; OUTPUTS:
;   nestmap  - healpix map in nested order
;
; COMMENTS:
;   
; REVISION HISTORY:
;   2003-Dec-04  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function map_ring2nest, ringmap

  npix = n_elements(ringmap)
  nside = round(sqrt(npix/12))
       
  dpf_nest2ring, nside, lindgen(12*nside*nside), ipring

  return, ringmap[ipring]
end
