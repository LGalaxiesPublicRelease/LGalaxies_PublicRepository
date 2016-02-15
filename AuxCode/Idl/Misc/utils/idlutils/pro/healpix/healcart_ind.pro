;+
; NAME:
;   healcart_ind
;
; PURPOSE:
;   Generate index list for healcart projection
;
; CALLING SEQUENCE:
;   ind = healcart_ind(data, /nest)
;
; INPUTS:
;   data  - array to determine size (not used for anything else)
;            if single element, interpret as Nside. 
;
; KEYWORDS:
;   nest  - return index list for nested ordering
;
; OUTPUTS:
;   ind   - index array to transform healpix to quasi-cartesian
;
; EXAMPLES:
;   ind = healcart_ind(kband, /nest)
;   atv, kband[ind]
;
; COMMENTS:
;   This routine returns indices for the "healcart" projection, 
;    which preserves the ring -> row mapping but stretches in
;    longitude to make a (nearly) Cartesian projection. 
;   This routine should only be used for examining healpix maps -- any
;    quantitative computations should only be done on the healpix sphere. 
;
;   Note: there are 4N-1 rows, not 4N
;
; REVISION HISTORY:
;   2003-Feb-12  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function healcart_ind, data, nest=nest

  if n_elements(data) gt 1 then begin 
     npix  = n_elements(data) 
     nside = round(sqrt(npix/12))
  endif else begin 
     if NOT keyword_set(data) then begin 
        print, 'ind=healcart_ind(data, nest=nest)'
        return, -1
     endif 
     nside = long(data)
     npix  = 12*nside*nside
  endelse 

  nx = 8*nside
  x = lindgen(nx)
  ind = lonarr(nx, 4*nside-1, /nozero)
  for i=0L, nside-2 do ind[*, i] = 2*i*(i+1)+4*(i+1)*x/nx

  xoffs = 2L*nside*(nside-1)
  for i=0L, 2*nside do ind[*, i+nside-1] = $
    shift(xoffs+4*nside*i+x/2, -(i mod 2))

  for i=1L, nside-1 do ind[*, 4*nside-1-i] = npix-2*i*(i+1)+4*i*x/nx

  ind = rotate(ind, 2)
  ind = shift(ind, 4*nside, 0)

  if keyword_set(nest) then begin 
     dpf_nest2ring, nside, lindgen(12*nside*nside), ipring
     ipnest = lonarr(nx, 4*nside-1, /nozero)
     ipnest[ipring] = lindgen(12*nside*nside)
     ind = ipnest[ind]
  endif 

  return, ind
end


