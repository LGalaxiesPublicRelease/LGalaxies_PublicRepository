;+
; NAME:
;   psf_findstars
;
; PURPOSE:
;   Find stars suitable for using in PSF estimation
;
; CALLING SEQUENCE:
;   psf_findstars, image, ivar, npad, clean, xstar, ystar, nsigma=, $
;                  satmask=, badpixels=
;
; INPUTS:
;   image      - image to locate stars in
;   ivar       - inverse variance image -- must be correctly
;                calibrated, or some of the algorithms will NOT WORK!
;
; OPTIONAL INPUTS:
;   npad       - do not consider npad pixels around the edge.  DEFAULT 2
;   
; KEYWORDS:
;   nsigma     - only consider stars brighter than nsigma sigma above mean
;   satmask    - saturation mask (1=bad)
;                badpixels  - bad pixel mask (should mask out
;                             everything around bad columns, bleeds, etc.)
;
; OUTPUTS:
;   clean      - clean image, at least clean enough to use for PSF
;                estimation. 
;   {x,y}star  - integer pixel locations of stars (0-indexed)
; 
; RESTRICTIONS:
;   This routine simply will not work right if you do not pass the
;   correct ivar array, because it calls psf_reject_cr(). 
;   (If you don't know your gain, we can't help you.)
;
; EXAMPLES:
;   
; COMMENTS:
;   This routine does NOT return an exhaustive list of stars.  It is
;    somewhat selective about avoiding the bad pixel mask, etc.  Also,
;    it will return a few CRs, if they are roughly band-limited. 
;    But it will reject most of the CRs (unless they are on the bad
;    pixel mask!)
;
;   After you have the PSF, you should call a real CR rejection
;    routine to properly clean the image and a real object finder. 
;
; REVISION HISTORY:
;   2006-May-26  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
pro psf_findstars, image, ivar, npad, clean, xstar, ystar, $
       nsigma=nsigma, satmask=satmask, badpixels=badpixels

; -------- set defaults
  if NOT keyword_set(image)     then message, 'must set image'
  if NOT keyword_set(ivar)      then message, 'must set ivar'
  if NOT keyword_set(npad)      then npad = 2
  if NOT keyword_set(nsigma)    then nsigma = 20
  if NOT keyword_set(satmask)   then message, 'you really should set satmask'
  if NOT keyword_set(badpixels) then message, 'works better if you set badpixels'
; -------- statistics (squash image for speed)
  sz = size(image, /dim)
  row = lindgen((sz[1]-1)/16+1)*16
  squash = image[*, row]
  squashbad = badpixels[*, row]
  wgoodsquash = where(squashbad EQ 0, ngoodsquash)
  if ngoodsquash LT 100 then message, 'your whole image is bad!'
  djs_iterstat, squash[wgoodsquash], sigma=sigma, mean=mean, sigrej=3

; -------- high pixels
  im  = (image-mean)*(1B-badpixels)
  ind = where(im GT (nsigma*sigma), nhigh)

  if nhigh eq 0 then message, 'no high pixels!'

  ix = ind mod sz[0]
  iy = ind / sz[0]

  ixm = (ix-1) > 0
  ixp = (ix+1) < (sz[0]-1)

  iym = (iy-1) > 0
  iyp = (iy+1) < (sz[1]-1)

  maxnb = im[ixm, iy] > im[ixp, iy] > im[ix, iym] > im[ix, iyp] > $
    im[ixm, iyp] > im[ixp, iyp] > im[ixp, iym] > im[ixm, iym]

  bd = (ivar EQ 0) OR badpixels
  badnb = bd[ixm, iy] OR bd[ixp, iy] OR bd[ix, iym] OR bd[ix, iyp] OR $
    bd[ixm, iyp] OR bd[ixp, iyp] OR bd[ixp, iym] OR bd[ixm, iym]

; -------- demand that the peak is larger than the neighbor
  peak = (im[ind] GT (maxnb*1.00001)) AND (badnb EQ 0)
  peak = peak AND ((satmask OR badpixels)[ind] EQ 0)

; -------- and not near the edge of the image...
  peak = peak AND ((iy GT npad) AND (iy LT (sz[1]-npad-1)) AND $
    (ix GT npad) AND (ix LT (sz[0]-npad-1)))

; -------- real peaks (mostly) - still need to worry about
;          diff. spikes, badcols
  w = where(peak, npeak)

; -------- we compute these ratios and then do not use them -- but
;           might want to in the future. 
; 02 12 22
; 01 11 21
; 00 10 20 

  im00 = im[ixm[w], iym[w]]
  im10 = im[ix[w],  iym[w]]
  im20 = im[ixp[w], iym[w]]

  im01 = im[ixm[w], iy[w]]
  im11 = im[ix[w],  iy[w]]
  im21 = im[ixp[w], iy[w]]

  im02 = im[ixm[w], iyp[w]]
  im12 = im[ix[w],  iyp[w]]
  im22 = im[ixp[w], iyp[w]]

  back1 = (im01+im21)/2.
  back2 = (im10+im12)/2.
  back = (back1+back2)/2.
  diag = (im00+im02+im20+im22)/4.

; -------- compute psfvals
  psfvals = [median(back/im11), median(diag/im11)]

  psfvals0 = [.46, .21]         ; assuming 1.88 pixel FWHM or bigger
                                ; if the PSF is smaller, you are
                                ; BADLY sampled!
  
; -------- fall back to marginal sampling psf if something has gone
;          wrong!
;  if (psfvals[0] LT psfvals0[0]) OR (psfvals[1] LT psfvals0[1]) then begin
  splog, 'median psfvals: ', psfvals
  psfvals=psfvals0
  splog, 'Using psfvals = ', psfvals

  rat = (((back) > 0) )/im11
  rat2 = ((diag > 0))/im11

  wstar = where(rat GT psfvals[0] and rat2 GT psfvals[1], nstar)

; -------- Find ALL CRs with psf_reject_cr
  cr = psf_reject_cr(image-mean, ivar, psfvals, satmask=satmask)
  clean = djs_maskinterp(image-mean, cr, iaxis=0, /const)

; -------- return only the good stars
  xstar = ix[w[wstar]]
  ystar = iy[w[wstar]]

; -------- report number of stars used
  splog, nstar, ' stars found'

  return
end
