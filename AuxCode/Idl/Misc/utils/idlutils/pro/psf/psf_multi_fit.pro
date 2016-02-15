;+
; NAME:
;   psf_multi_fit
;
; PURPOSE:
;   Fit multiple PSFs to a postage stamp, given a preliminary PSF.
;
; CALLING SEQUENCE:
;   psf_multi_fit, stamps, stampivar, psfs, par, sub, faint, nfaint=nfaint
;
; INPUTS:
;   
; OPTIONAL INPUTS:
;   
; KEYWORDS:
;   
; OUTPUTS:
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
;   2006-May-25   Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function faint_struc, N

  nmax = 5
  str = {npsf:0, $
         px:fltarr(nmax), $
         py:fltarr(nmax), $
         peak:fltarr(nmax)}
  
  if keyword_set(N) then begin
     str = replicate(str, N)
  endif 

  return, str
end



pro psf_multi_fit, stamps, stampivar, psfs, par, sub, faint, nfaint=nfaint

  sub = stamps                        ; copy of array for output
  nstamp = (size(stamps, /dimen))[2]
  npix = (size(stamps, /dimen))[0]
  nsigma = 5
  faint = faint_struc(nstamp)
  if NOT keyword_set(nfaint) then nfaint = 3
  if nfaint GT n_elements(faint[0].peak) then message, 'nfaint too large'
  rad = 3

  xbox = djs_laxisnum([npix, npix], iaxis = 0)
  ybox = djs_laxisnum([npix, npix], iaxis = 1)

  mask = bytarr(npix, npix)
  npad = par.boxrad-par.fitrad
  mask[npad:npix-npad-1, npad:npix-npad-1] = 1B

  xcen = (npix-1)/2
  ycen = (npix-1)/2
  mask = mask AND (((xbox-xcen)^2+(ybox-ycen)^2) GT rad^2)

; -------- loop over stamps
  for i=0L, nstamp-1 do begin 

     fmask = mask

     psf = psfs[*, *, i]
     im = stamps[*, *, i]-psf
     iv = stampivar[*, *, i]
     status = 1B

     for j=0, nfaint-1 do begin 
; -------- look for another peak
        maxval = max(im*fmask, maxind)
        if (maxval*sqrt(iv[maxind]) GT nsigma) and (status EQ 1B) then begin 

           ix = maxind mod npix
           iy = maxind  /  npix
        
; -------- get sub-pixel center

           junk = psf_stamp_center_iter(im, 1, dx=dx0, dy=dy0, center=[ix, iy], maxiter=2, status=status)

           if ((abs(dx0) > abs(dy0)) GE 2) or (status NE 1B) then begin 
              print, 'Centering error on stamp', i, dx0, dy0
           endif else begin 
              px = ix+dx0
              py = iy+dy0

; -------- subtract shifted PSF form im
              spsf = sshift2d(psf, [px, py]-(npix-1)/2)
              norm = maxval/max(spsf)
              im = im-spsf*norm
              
; -------- store results
              faint[i].npsf = faint[i].npsf+1
              faint[i].px[j] = px
              faint[i].py[j] = py
              faint[i].peak = norm

              fmask = fmask AND ((xbox-px)^2+ $
                                 (ybox-py)^2) GT rad^2

           endelse
        endif
     endfor

     sub[*, *, i] = im
  endfor
  return
end
