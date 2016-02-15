;+
; NAME:
;   psf_polyfit
;
; PURPOSE:
;   Fit polynomial as a function of image (x,y) to each pixel in stamp
;
; CALLING SEQUENCE:
;   cf = psf_polyfit(stack, ivar, x, y, par, ndeg=ndeg, reject=reject, $
;            cond=cond, scale=scale)
; INPUTS:
;   stack    - postage stamps to fit (npix, npix, nstar)
;   ivar     - ivar for same
;   {x,y}    - star positions ([0..2048,0..1361] for SDSS)
;   par      - par structure
;   ndeg     - degree of fit
;
; OPTIONAL INPUTS:
;   scale    - range of (x,y) (2 element array) to rescale (x,y) to [-1,1]
;
; KEYWORDS:
;   reject   - set to reject stars that don't fit well and refit
;
; OUTPUTS:
;   ndeg     - degree of fit ACTUALLY DONE.
;
; OPTIONAL OUTPUTS:
;   cond     - worst condition number encountered in fit
;
; RESTRICTIONS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   Outliers are rejected on a star by star basis, not just pixel by
;    pixel. 
;   If there are not enough stars for the requested degree, falls back
;    to a more appropriate number. 
;   Calling routine should inspect Condition number and reduce degree
;    further if necessary.
;
;   We should have a noise floor in the ivar!!!
;   
; REVISION HISTORY:
;   2006-May-27   Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function psf_polyfit, stack, ivar, x, y, par, ndeg=ndeg, reject=reject, $
            cond=cond, scale=scale

; -------- check inputs
  psfsize = (size(stack, /dimen))[0]
  if psfsize NE (2*par.boxrad+1) then message, 'parameter mismatch'
  if NOT keyword_set(scale) then scale = double([2048, 1361])
  if n_elements(scale) NE 2 then message, 'scale should have 2 elements'
  if size(stack, /n_dimen) EQ 2 then begin
     cf = stack
     ndeg = 0
     cond = cf*0+1.0
     return, cf
  endif

  x0 = par.boxrad-par.fitrad
  x1 = par.boxrad+par.fitrad

  nstamp = (size(stack, /dimen))[2]
  ncoeff = (ndeg+1)*(ndeg+2)/2
  if n_elements(ndeg) EQ 0 then ndeg = 1

; -------- check that ndeg isn't too high
  if ncoeff GE nstamp then begin 
     ndeg = floor((sqrt(1+8*(nstamp-1))-3)/2)
     ncoeff = (ndeg+1)*(ndeg+2)/2
     splog, 'Falling back to ndeg = ', ndeg
  endif 

  cf   = fltarr(psfsize, psfsize, ncoeff)
  cond = fltarr(psfsize, psfsize)

; -------- construct A matrix
  A = dblarr(ncoeff, nstamp)
  col = 0

; -------- map (x,y) to (-1.1) so condition number makes sense
  xd = double(x)/scale[0]*2-1
  yd = double(y)/scale[1]*2-1
  for ord=0L, ndeg do begin
     for ypow=0, ord do begin 
        xpow = (ord-ypow)
        A[col, *] = xd^xpow * yd^ypow
        col = col+1
     endfor
  endfor

; -------- do fit
;  W = dblarr(nstamp)+1
  for i=x0, x1 do begin 
     for j=x0, x1 do begin 
        W = reform(ivar[i, j, *])        ; should have noise floor here!

        W = W < (min(W)*100)
        data = reform(stack[i, j, *])
        hogg_iter_linfit, A, data, W, coeff, nsigma=5, /median, $
          /truesigma, condition=condition

        cf[i, j, *] = coeff
        cond[i, j] = condition
     endfor
  endfor

; -------- reject
  if keyword_set(reject) then begin 
     nsigma = 4
     model = psf_eval(x, y, cf, par.cenrad)
     bad = (model-stack)^2*ivar GT nsigma^2
     nbad = total(total(bad[x0:x1, x0:x1, *], 1), 1)
     good = where(nbad LE 5, ngood)
     if ngood LT 1 then message, 'Bad news - we rejected all the stars!'

     splog, 'keeping ', ngood, ' stars.'
     cf = psf_polyfit(stack[*, *, good], ivar[*, *, good], x[good], y[good], $
                  par, ndeg=ndeg, cond=cond, scale=scale)

  endif

  return, cf
end
