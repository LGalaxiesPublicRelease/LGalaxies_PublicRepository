;+
; NAME:
;   psf_eval
;
; PURPOSE:
;   Evaluate PSF at (x,y), possibly sinc-shifted to center (dx,dy)
;
; CALLING SEQUENCE:
;   psfs = psf_eval(x, y, coeff, cenrad, dx=dx, dy=dy, scale=scale)
;
; INPUTS:
;   x,y     - image coordinates at which to evaluate PSF
;   coeff   - coefficients from psf_polyfit()
;   cenrad  - radius used to center and normalize.  From par struct.
;
; OPTIONAL INPUTS:
;   dx,dy   - subpixel shift for sinc-shift
;   scale   - range of (x,y) (2 element array) to rescale (x,y) to [-1,1]
;
; OUTPUTS:
;   psfs    - (npix, npix, nstar) array of psfs
;
; COMMENTS:
;   see psf_polyfit()
;
; REVISION HISTORY:
;   2006-May-27   Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function psf_eval, x, y, coeff, cenrad, dx=dx, dy=dy, scale=scale

  psfsizex = (size(coeff, /dimen))[0]
  psfsizey = (size(coeff, /dimen))[1]
  if size(coeff, /n_dimen) EQ 3 then begin
     ncoeff = (size(coeff, /dimen))[2]
  endif else ncoeff = 1
  npsf = n_elements(x) 
  ndeg = long(sqrt(ncoeff*2))-1
  if (ndeg+1)*(ndeg+2)/2 NE ncoeff then stop
  if NOT keyword_set(scale) then scale = double([2048, 1361])
  if n_elements(scale) NE 2 then message, 'scale should have 2 elements'

  
; -------- compute A matrix
  A = dblarr(ncoeff, npsf)
  col = 0
  xd = double(x)/scale[0]*2-1
  yd = double(y)/scale[1]*2-1
  for ord=0L, ndeg do begin
     for ypow=0, ord do begin 
        xpow = (ord-ypow)
        A[col, *] = xd^xpow * yd^ypow
        col = col+1
     endfor
  endfor

  psfs = fltarr(psfsizex, psfsizey, npsf)
  for i=0L, psfsizex-1 do begin 
     for j=0L, psfsizey-1 do begin 
        psfs[i, j, *] = reform(coeff[i, j, *])#A
     endfor
  endfor

  if keyword_set(dx) AND keyword_set(dy) then begin 
     for i=0L, npsf-1 do begin 
        psfs[*, *, i] = sshift2d(psfs[*, *, i], [dx[i], dy[i]])
        psfs[*, *, i] = psfs[*, *, i]/psf_norm(psfs[*, *, i], cenrad)
     endfor
  endif

  return, psfs
end
