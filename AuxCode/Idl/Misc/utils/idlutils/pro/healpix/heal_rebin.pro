;+
; NAME:
;   heal_rebin
;
; PURPOSE:
;   rebin ring ordered healpix maps to different Nside
;
; CALLING SEQUENCE:
;   newmap = heal_rebin(oldmap, newnside)
;
; INPUTS:
;   oldmap   - original map, ring order
;
; OUTPUTS:
;   newmap   - new rebinned map, ring order
;
; EXAMPLES:
;   
; COMMENTS:
;   Only call with ring-ordered maps.
;   Use IDL rebin() if your map is already nested order. 
;
; REVISION HISTORY:
;   2003-Dec-04  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function heal_rebin, oldmap, newnside

  newnpix = 12L*newnside*newnside
  oldnpix = n_elements(oldmap) 
  if oldnpix EQ newnpix then return, oldmap

  oldnest = map_ring2nest(oldmap)
  if newnpix LT oldnpix then newnest = rebin(oldnest, newnpix)
  if newnpix GT oldnpix then newnest = rebin(oldnest, newnpix, /sample)

  newmap = map_nest2ring(newnest)

  return, newmap
end
