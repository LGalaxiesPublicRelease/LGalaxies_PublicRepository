;------------------------------------------------------------------------------
;+
; NAME:
;   djs_phot
;
; PURPOSE:
;   Driver for aperture photometry with the option of re-centering and
;   sky-subtraction.
;
; CALLING SEQUENCE:
;   flux = djs_phot( xcen, ycen, objrad, skyrad, image, [invvar, $
;    calg=, cbox=, cmaxiter=, cmaxshift=, $
;    fwhm=, fixfw=, ceps=, $
;    salg=, srejalg=, smaxiter=, $
;    lorej=, hirej=, $
;    flerr=, skyval=, skyrms=, skyerr=, peakval=, /quick, /exact ] )
;
; INPUTS:
;   xcen:       X center(s)
;   ycen:       Y center(s)
;   objrad:     Radius for aperture on object, or a vector of such radii.
;   skyrad:     A 2-element array with two radii to define an annulus,
;               or a scalar to define a circular aperture, or undefined
;               for no sky calculation
;   image:      FITS data array, as read from readfits().
;
; OPTIONAL INPUTS:
;   invvar:     Inverse variance image, for computing the errors FLERR
;   ----------- FOR CENTERING ALGORITHM
;   calg:       Centering algorithm.  Choose from iweight, gauss1, gauss2, none.
;                 iweight = intensity-weighted center, computed independently
;                           in both X and Y
;                 gauss1  = gaussian fit, including a constant term, computed
;                           independently in both X and Y
;                 gauss2  = not implemented
;                 none    = no centering
;               Default to iweight.
;   cbox:       Centering box width.  Default to 7.
;   cmaxiter:   Maximum number of iterations for centering algorithm.
;               Default to 10.
;   cmaxshift:  Maximum center shift.  If this shift is exceeded in either
;               X or Y, then return the center position passed in XCEN,YCEN.
;               A value of 0 imposes no maximum shift.  Default to 3.
;   fwhm:       FWHM for gauss1 or gauss2 centering algorithms.  If this is
;               a scalar, then the same value is used for both X and Y.
;               If this is a vector, then the first two elements are used
;               for X and Y respectively.
;   fixfw:      If set and nonzero, then fix the FWHM for gauss1 or gauss2 fits.
;   ceps:       Stop iterating when relative shifts in both X and Y are less
;               than this value.  Default to 0.
;
;   ----------- FOR SKY FITTING ALGORITHM
;   salg:       Sky fitting algorithm.  Choose from mean, median, mode, none.
;               Default to "mean" if SKYRAD is set, or "none" otherwise.
;   srejalg:    Rejection algorithm.  Choose from none, sigclip, pclip.
;                 sigclip = sigma clipping; reject outliers that are
;                           more than lorej*sigma below skyval or hirej*sigma
;                           above skyval
;                 pclip   = percentile clipping; reject the lowest lorej
;                           fraction and the highest hirej fraction of points
;                 none    = no rejection
;               Default to sigclip
;   smaxiter:   Maximum number of srejalg iterations.  Default to 10.
;   lorej:      If srejalg="sigclip", then the number of standard deviations
;               below skyval to clip (default to 3.0).
;               If srejalg="pclip", then fraction of low pixels to clip
;               (default to 0.05).
;   hirej:      If srejalg="sigclip", then the number of standard deviations
;               above skyval to clip (default to 3.0).
;               If srejalg="pclip", then fraction of high pixels to clip
;               (default to 0.05).
; 
; KEYWORDS:
;   exact       Use slow photo-based photfrac algorithm (Not thoroughly tested)
;               Requires image to be centered such that xcen and ycen
;               are integer values. If set, does not recalculate
;               center.
;   quick       Use faster photfrac algorithm (Not thoroughly tested)
;
; OUTPUTS:
;   flux:       Total flux within each circular aperture defined by objrad,
;               minus the sky contribution within each aperture [NOBJ,NRAD].
;   xcen:       Re-centered X position (modified).
;   ycen:       Re-centered X position (modified).
;
; OPTIONAL OUTPUTS:
;   flerr:      Flux error from the sky background uncertainty within
;               the object aperture(s) [NOBJ,NRAD]
;   skyval:     Sky value in counts per unit area [NOBJ].
;   skyrms:     RMS of sky pixels about skyval, again in counts per unit area.
;               This assumes that each unit area is independent.  The RMS
;               is computed for only the values that remain after all the
;               rejection iterations [NOBJ].
;   skyerr:     The error in skyval assuming that each pixel is independent.
;               This is skyrms divided by the square root of the number of
;               pixels [NOBJ].
;   peakval:    Peak pixel value (before sky-subtraction)
;
; COMMENTS:
;   Sub-pixel sampling of the circular apertures is handled exactly.
;   If /exact keyword is set, input xcen, ycen must be integers or 
;     the code bombs. See exact_photfrac.pro for more details, but 
;     basically exact_photfrac is much simpler to implement if the 
;     object is already centered, which doesn't cost you precision.
;   For similar reasons, if /exact is set, djs_phot will not try to
;     recentroid your object.
; PROCEDURES CALLED:
;   djs_photcen
;   djs_photfrac
;   djs_photsky()
;
; REVISION HISTORY:
;   28-Nov-1996  Written by D. Schlegel, Durham.
;   01-Jun-2000  Major revisions: change XYCEN to XCEN,YCEN; remove use
;                of FITS headers; make IDL 5 compliant (DJS).
;   02-Nov-2000  objrad, skyrad recast as floats (D. Finkbeiner)
;                  If they are ints, 1% errors may arise. 
;-
;------------------------------------------------------------------------------
function djs_phot, xcen, ycen, objrad, skyrad, image, invvar, $
 calg=calg, cbox=cbox, cmaxiter=cmaxiter, cmaxshift=cmaxshift, $
 fwhm=fwhm, fixfw=fixfw, ceps=ceps, $
 salg=salg, srejalg=srejalg, smaxiter=smaxiter, $
 lorej=lorej, hirej=hirej, $
 flerr=flerr, skyval=skyval, skyrms=skyrms, skyerr=skyerr, peakval=peakval, $
 quick=quick, exact=exact

   ; Need 5 parameters
   if (N_params() LT 5) then begin
      doc_library, 'djs_phot'
      return, -1
   endif
   nobj = n_elements(xcen)
   if (n_elements(ycen) NE nobj) then $
    message, 'XCEN and YCEN must have same number of elements'

   if(keyword_set(exact)) then begin
       ini=where(long(xcen) ne xcen OR long(ycen) ne ycen, nni)
       if(nni gt 0) then begin
           message, 'xcen and ycen MUST be integers if /exact keyword is set'
       endif
   endif

   dims = size(image, /dimens)
   xdimen = dims[0]
   ydimen = dims[1]
   nrad = n_elements(objrad)

   ; Allocate memory for return values
   flux = fltarr(nobj,nrad)
   if (arg_present(flerr)) then flerr = fltarr(nobj,nrad)
   skyval = fltarr(nobj)
   skyrms = fltarr(nobj)
   skyerr = fltarr(nobj)
   if (arg_present(peakval)) then peakval = fltarr(nobj,nrad)

   ;----- LOOP THROUGH EACH OBJECT -----
   for iobj=0L, nobj-1 do begin

      ; Center the object
      xcen1 = xcen[iobj]
      ycen1 = ycen[iobj]
      if(NOT keyword_set(exact)) then begin
          djs_photcen, xcen1, ycen1, image, $
            calg=calg, cbox=cbox, cmaxiter=cmaxiter, cmaxshift=cmaxshift 
      endif else begin
          splog, 'Not recentering -- /exact requires correct input center'
          xcen1=xcen[iobj]
          ycen1=ycen[iobj]
      endelse

      ; Find the sky value
      ; Add 0.0 to the skyrad to convert to floating-point
      skyval[iobj] = djs_photsky( xcen1, ycen1, skyrad+0.0, image, $
       salg=salg, srejalg=srejalg, smaxiter=smaxiter, $
       lorej=lorej, hirej=hirej, $
       skyrms=tmprms, skyerr=tmperr, quick=quick, exact=exact)
      skyrms[iobj] = tmprms
      skyerr[iobj] = tmperr

      ; Find total counts in object aperture(s)
      for irad=0L, nrad-1 do begin
         onerad = objrad[irad] + 0.0 ; Convert to floating-point
         if (onerad GT 0) then begin
            if keyword_set(quick) then begin 
               quick_photfrac, xcen1, ycen1, onerad, $
                 xdimen=xdimen, ydimen=ydimen, $
                 pixnum=pixnum, fracs=fracs
           endif else if keyword_set(exact) then begin 
               exact_photfrac, xcen1, ycen1, onerad, $
                 xdimen=xdimen, ydimen=ydimen, $
                 pixnum=pixnum, fracs=fracs
            endif else begin 
               djs_photfrac, xcen1, ycen1, onerad, $
                 xdimen=xdimen, ydimen=ydimen, $
                 pixnum=pixnum, fracs=fracs
            endelse
         endif else begin
             pixnum = -1L
         endelse
         if (pixnum[0] EQ -1) then begin
            flux[iobj,irad] = 0.0  ; No object pixels
            if (arg_present(flerr)) then flerr[iobj,irad] = 0.0
         endif else begin
            area = total(fracs)
            flux[iobj,irad] = total( image[pixnum] * fracs ) $
             - (area * skyval[iobj])
            if (arg_present(flerr)) then flerr[iobj,irad] = area * skyerr[iobj]
            if (arg_present(peakval)) then $
             peakval[iobj,irad] = max(image[pixnum])
         endelse

         ; If INVVAR was set, then measure those pixel errors and add
         ; them in quadrature to the sky-subtraction errors.
         if (keyword_set(invvar) AND arg_present(flerr)) then begin
            if (pixnum[0] EQ -1) then sigma2 = 0 $
             else if (min(invvar[pixnum]) LE 0) then sigma2 = 0 $
              else sigma2 = total( fracs / invvar[pixnum] )
            if (sigma2 LE 0) then $
             flerr[iobj,irad] = -1 $
            else $
             flerr[iobj,irad] = sqrt(flerr[iobj,irad]^2 + sigma2)
         endif
      endfor

      ; Overwrite the center positions with the newly computed ones
      xcen[iobj] = xcen1
      ycen[iobj] = ycen1

   endfor

   return, flux
end
;------------------------------------------------------------------------------
