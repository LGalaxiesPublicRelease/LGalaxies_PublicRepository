;+
; NAME:
;   psf_reject_cr
;
; PURPOSE:
;   Find pixels contaminated by cosmic rays (CRs) using PSF criteria
;
; CALLING SEQUENCE:
;   crmask = psf_reject_cr(image, ivar, [ psfvals, ] satmask=, $
;    [ pstr=, nsigma=, cfac=, niter=, c2fudge= ] )
;
; INPUTS:
;   image     - image to test
;   ivar      - inverse variance of image
;   psfvals   - values of psf at 1 pix and sqrt(2) pixels from center
;               [2-element array]; must be set if PSTR is not
;   satmask   - saturated pixel mask (1=saturated)
;
; OPTIONAL KEYORDS:
;   pstr      - ???
;   nsigma    - minimum value (in sigma) for condition 2 [default 6]
;   cfac      - consistency factor [in sigma] for condition 3 [default 3]
;   niter     - max number of iterations [default 6]
;   c2fudge   - fudge factor to multiply psfvals by on first pass [default 0.8] 
;
; OUTPUTS:
;   crmask    - byte array, same size as image (1=CR)
;
; EXAMPLES:
;
; COMMENTS:
;   Algorithms designed by R Lupton and J. Gunn, implemented in C by Lupton,
;    re-implemented by M. Blanton as reject_cr. 
;    Now completely rewritten by D. Finkbeiner as psf_reject_cr.
;    see Lupton's CR.c in photo product for more. 
;
;   Do NOT need to call with zeroed image.
;   Lupton's screed is at: http://www.astro.princeton.edu/~rhl/photomisc/
;
; REVISION HISTORY:
;   2005-Mar-09  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function psf_reject_cr, image, ivar, psfvals, satmask=satmask, $
 nsigma=nsigma, cfac=cfac, niter=niter, c2fudge=c2fudge, pstr=pstr

; -------- check inputs
  if NOT keyword_set(image) then message, 'must set image'
  if NOT keyword_set(ivar) then message, 'must set ivar'
  if NOT (keyword_set(psfvals) OR keyword_set(pstr)) then $
   message, 'Must set PSFVALS or PSTR'
  if NOT keyword_set(nsigma) then nsigma = 6.0
  if NOT keyword_set(cfac) then cfac = 3.0
  if NOT keyword_set(niter)  then niter = 6
  if NOT keyword_set(satmask) then message, 'need to set satmask'
  if NOT keyword_set(c2fudge) then c2fudge = 0.8

; -------- make copy of image to work on
  im = image

  sz = size(im, /dimens)
  gd = bytarr(sz[0], sz[1])+1B
  cr = bytarr(sz[0], sz[1])

; -------- sample some rows and take minimum sig and sky
  nrow = 100
  drow = sz[1]/nrow
  yrow = lindgen(nrow)*drow
  sky = 1e+38
  sig = 1e+38
  for i=0, nrow-1 do begin 
     djs_iterstat, im[*, yrow[i]], mean=sky0, sigma=sig0, sigrej=3
     sky = sky < sky0
     sig = sig < sig0
  endfor

  im = im-sky
  msig = im GT (nsigma*sig)

; -------- suspects (indexes of im) are above threshold on the first
;          pass, and neighbors of CRs on subsequent passes
  suspects = where(msig, nsuspect)
  if nsuspect EQ 0 then return, cr  ; cannot find anything

  for iiter=0L, niter-1 do begin 
                                ; cr_find investigates suspects,
                                ; returns mpeak of same lenght (1=CR, 0=OK)

; -------- if pstr set then get psfvals for each "suspect"
     if keyword_set(pstr) then begin 
        x0 = suspects mod sz[0]
        y0 = suspects / sz[0]
        cf = pstr.coeff
        cen = pstr.boxrad
        psfcore = psf_eval(x0, y0, cf[cen-1:cen+1, cen-1:cen+1, *])
        peak = reform(psfcore[1, 1, *])
        psfval1 = reform((psfcore[0, 1, *]+psfcore[2, 1, *])/2)/peak
        psfval2 = reform((psfcore[1, 2, *]+psfcore[1, 0, *])/2)/peak
        psfval3 = reform((psfcore[0, 2, *]+psfcore[2, 0, *])/2)/peak
        psfval4 = reform((psfcore[0, 0, *]+psfcore[2, 2, *])/2)/peak
        psfvals_use = transpose([[(psfval1+psfval2)/2.], [(psfval3+psfval4)/2.]])
     endif else psfvals_use = psfvals


     c3fac = iiter EQ 0 ? cfac : 0.0
     thisc2 = iiter EQ 0 ? c2fudge[0] : 1.0
     mpeak = psf_reject_cr_single(im, gd, ivar, satmask, nsigma, $
                                  c3fac, psfvals_use*thisc2, $
                                  suspects, neighbor=neighbor)
     wrej  = where(mpeak, nrej)

; -------- escape if we are done     
     if nrej EQ 0 then return, cr
  
; -------- otherwise examine neighbors
     cr[suspects[wrej]] = 1B
     gd[suspects[wrej]] = 0B
     suspects = neighbor

  endfor

  return, cr
end
