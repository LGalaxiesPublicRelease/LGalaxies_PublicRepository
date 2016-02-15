;+
; NAME:
;   find_sb
; PURPOSE:
;   Find objects in a 2-d image with a gaussian filter
; CALLING SEQUENCE
;   find_sb, sub, wt, x=x, y=y, flux=flux, $
;       [sigma=, hpix=, hmin=]
; INPUTS:
;   sub         - skysubtracted image 
;   wt          - inverse variance weight (0 or 1 for CR mask is OK)
; OPTIONAL INPUTS:
;   sigma       - sigma of gaussian filter in pixels  (1.0 is default)
;   hpix        - half-pixel width of convolution filter  (2 is default)
;   hmin        - minimum flux threshold  
;   curvr       - the maximum allowed log ratio of curvature along the 
;                     major/minor axes  (basically checking for roundness).
;                   default is 2.
; KEYWORDS:
; OUTPUTS:
;   x          - column pixel positions of flux peak (sorted be decreasing flux)
;   y          - row pixel position                  "
;   flux       - gaussian filtered flux estimate
;   pa_degrees - a guess at the PA of the object from the x-axis
;   ab         - a guess at the minor/major axis ratio
; COMMENTS:
; DEPENDENCIES:
;   idlutils
; BUGS:
;   No checks on neighboring peaks
;   Not tested yet with real inverse variance weighting
;   Not sure that I have PA calculated correctly
; REVISION HISTORY:
;   2005-11-30  Written - Burles
;-
pro find_sb, sub, wt, x=x, y=y, flux=flux, sigma=sigma, $
             pa_degrees=pa_degrees, ab=ab, fwhm=fwhm, $
             hpix=hpix, hmin=hmin, curvr=curvr

 
  if NOT keyword_set(sigma) then sigma = 1.0
  if NOT keyword_set(hpix ) then hpix = 2L
  if NOT keyword_set(curvr ) then curvr = 2.0

  npix = 2*hpix + 1L
  k = gauss_kernel(sigma, hpix=2) 
  k2 = k # k

  denom_opt = convol(wt, k2^2)
  denom_ap  = convol(wt, k2)
  numer = convol(sub*wt, k2) 
  f_opt = numer / (denom_opt + (denom_opt EQ 0))
  f     = numer / (denom_ap + (denom_ap EQ 0)) 
  frac = convol(1.0*(wt GT 0),k2)

;  mask = smooth(1.0*(wt LE 0),3) EQ 0   ; grow the wt EQ 0 mask
;  h = f*mask
;  if max(h) LE 0 then return

  if NOT keyword_set(hmin) then $
     hmin = 20. * median(abs(f - median(f)))

  index = where(f * frac GE hmin, nfound)
  
  n_x = (size(sub))[1]
  x = ((lindgen(npix^2) mod npix) - hpix)
  y = ((lindgen(npix^2) / npix) - hpix)
  off = long(x + y*n_x)
  r = sqrt(x^2 + y^2)

  lim = (2.0*sigma) > 2
  good = where(r LE lim AND r GT 0, npixels)
  offset = (good mod npix) - hpix + ((good /npix) - hpix) * n_x

  for i=0, n_elements(offset)-1L do begin & $
     stars = where( f[index] GE  f[index+offset[i]], nfound) & $
     if nfound LE 0 then break & $
     index = index[stars] & $
  endfor

  print, 'Checking ', nfound, ' peaks'

  spots = index ## replicate(1,npix^2) + off # replicate(1,nfound)
  im = f[spots]
  wt_im = wt[spots]

  para = fit_para2d(x,y,im,wt_im, model=model)
  pt = transpose(para)
  sq = sqrt((pt[*,3] -pt[*,5])^2 + pt[*,4]^2)
  l1 = pt[*,3]  + pt[*,5] + sq
  l2 = pt[*,3]  + pt[*,5] - sq
  det = l1*l2
  r1 = l1/(pt[*,0] + (pt[*,0] EQ 0))
  r2 = l2/(pt[*,0] + (pt[*,0] EQ 0))
  b = -1.0*para[1:2,*]
  xs = total(b * [[2*para[5,*], -para[4,*]]],1) / (det + (det EQ 0))
  ys = total(b * [[-para[4,*], 2*para[3,*]]],1) / (det + (det EQ 0))
  pa = atan(pt[*,4], pt[*,5]-pt[*,3])/2.      ; in radians

  keep = where(pt[*,0] GT hmin AND $
               (r1 GT -0.4)  AND $
               (r2 GT -0.4)  AND $
               l1*l2 GT 0.01*hmin^2 AND $
               abs(r1 -r2) LT 0.1 AND $
               alog(abs((l2+(l2 EQ 0))/(l1 + (l1 EQ 0)))) LT curvr AND $
               sqrt(xs^2 + ys^2) LT 0.71, nkeep)


  print, 'Keeping ', nkeep, ' peaks'

  ik = index[keep]
  x = (ik mod n_x) + xs
  y = (ik / n_x) + ys
  flux = f[ik] * ( 1.0 - (para[3,keep] * xs[keep]^2 + $
                          para[4,keep] * xs[keep] * ys[keep] + $
                          para[5,keep] * ys[keep]^2 )/para[0,keep])

  flux = flux / max(k2) *2. * !Pi * sigma^2
  s_f = reverse(sort(flux))

  x = x[s_f]
  y = y[s_f]
  flux = flux[s_f] 
  pa_degrees = 90. - 180*pa[keep[s_f]]/!Pi
  ab = (l1[keep]/l2[keep])[s_f]
  fwhm = 2.3548/2./sqrt(-r1[keep[s_f]])
  
  return
end
