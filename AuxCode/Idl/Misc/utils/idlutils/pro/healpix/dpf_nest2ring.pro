;+
; NAME:
;   dpf_nest2ring
;
; PURPOSE:
;   Convert healpix nested pixel number to ring pixel number
;
; CALLING SEQUENCE:
;   dpf_nest2ring, nside, ipnest, ipring
;
; INPUTS:
;   nside  - healpix nside
;   ipnest - nested pixel number
;
; OUTPUTS:
;   ipring - ring pixel number
;
; EXAMPLES:
;   
; COMMENTS:
;   Usage is same as Hivon's nest2ring, only this is 4 times as fast.
;
; REVISION HISTORY:
;   2003-Dec-04  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
pro dpf_nest2ring, nside, ipnest, ipring

  common dpf_nest2ring, dx, dy, nside_save

; -------- initialize parameters
  nside  = long(nside)
  nint   = fix(nside)
  npix   = 12L*nside*nside
  nlevel = round(alog(nside)/alog(2))
  if 2L^nlevel NE nside then message, 'bad nside value - must be power of 2.'

  maketable = 1B
  if keyword_set(nside_save) then begin 
     if nside_save eq nside then maketable = 0B
  endif 

  if maketable then begin 
; -------- make i0,i1 lookup
     p2 = 1L
     i0 = intarr(npix/12)
     i1 = intarr(npix/12)
     bits = fix(ipnest)

     isquare = lindgen(npix/12)
     for ilevel=0, nlevel-1 do begin
        thisbit = 2^ilevel
        i0 =  i0 OR (bits AND thisbit)
        bits = fix(ishft(isquare, -(1+ilevel)))
        i1 = i1 OR (bits AND thisbit)
     endfor

; -------- dx,dy defined for one square of the projection; need to be
;          indexed off of dind when used
     dy = i0+i1
     dx = temporary(i0)-temporary(i1)

; -------- cache dx, dy
     nside_save = nside
  endif 

  dind = ipnest AND (nside*nside-1)

; -------- square number, ring number
  sn      = fix(ishft(ipnest, -2*nlevel))
  ringnum = ishft(sn, -2)*nint-dy[dind]+(2*nint-2)

; -------- offsets along ring as a function of square number
  off = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7]

  npix_per_ring = [(indgen(nint)+1)*4, intarr(2*nint-1)+4*nint, $
                   (nint-indgen(nint))*4]
  nq = (npix_per_ring/4)[ringnum]

; -------- pixel offsets along each ring
  ringoff = ishft((dx[dind]+off[sn]*nq) AND (8*nint-1), -1)

  qind  = lindgen(nside)
  qindrev = nside-qind
  xoffs = 2L*nside*(nside+1)

; -------- ring0 is the lowest pixel number in each ring
  ring0 = [2*qind*(qind+1), xoffs+(4*nside)*lindgen(2*nside-1), $
           npix-2*qindrev*(qindrev+1)]

; -------- ring pixel index
  ipring = ring0[ringnum]+ringoff

  return
end
