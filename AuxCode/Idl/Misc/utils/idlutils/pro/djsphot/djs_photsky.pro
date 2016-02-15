;------------------------------------------------------------------------------
;+
; NAME:
;   djs_photsky
;
; PURPOSE:
;   Compute the sky value within an annulus.
;
;   At present, fractional pixels are not treated properly; this is because
;   we need a sort routine that carries index numbers, such as NR indexx().
;
; CALLING SEQUENCE:
;   skyval = djs_photsky( xcen, ycen, skyrad, image, $
;    [ salg=, srejalg=, smaxiter=, $
;    lorej=, hirej=, skyrms=, skyerr= ] )
;
; INPUTS:
;   xcen:       X center(s)
;   ycen:       Y center(s)
;   skyrad:     A 2-element array with two radii to define an annulus,
;               or a scalar to define a circular aperature, or undefined or 0
;               for no sky calculation
;   image:      FITS data array, as read from readfits().
;
; OPTIONAL INPUTS:
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
;   quick:      Set to use quick_photfrac (much faster)
;   exact:      Set to use exact_photfrac (slower)
;
; OUTPUTS:
;   skyval:     Sky value in counts per pixel.
;   skyrms:     RMS of sky pixels about skyval, again in counts per pixel.
;               This assumes that each pixel is independent.  The RMS
;               is computed for only the values that remain after all the
;               rejection iterations.
;   skyerr:     The error in skyval assuming that each pixel is independent.
;               This is skyrms divided by the square root of the number of
;               pixels.
;
; PROCEDURES CALLED:
;   djs_photfrac
;
; REVISION HISTORY:
;   28-Nov-1996  Written by D. Schlegel, Durham.
;   01-Jun-2000  Major revisions: change XYCEN to XCEN,YCEN; remove use
;                of FITS headers; make IDL 5 complient (DJS).
;-
;------------------------------------------------------------------------------
; INTERNAL SUPPORT PROCEDURES:
;
; djs_photsky_compute
; djs_photsky_reject
;------------------------------------------------------------------------------
function djs_photsky_compute, image, salg

   case salg of
   'mean': begin
      skyval = total(image) / N_elements(image)
      end

   'median': begin
      skyval = median(image)
      end

   'mode': begin
      skyval = 0.0
; NOT IMPLEMENTED!!!???
      end

   'none': begin
      skyval = 0.0
      end

   else: begin
      print, 'Error - Invalid sky algorithm'
      return, 0.0
      end
   endcase

   return, skyval
end
;------------------------------------------------------------------------------
function djs_photsky_reject, skyval, image, srejalg, lorej, hirej

   case srejalg of
   'sigclip': begin
      diff = image - skyval
      sigval = sqrt( total(diff*diff) / N_elements(diff) )
      pix = where(image GT (skyval-lorej*sigval) $
              AND image LT (skyval+hirej*sigval) )
      if (pix(0) EQ -1) then $
       newimg = image $
      else $
       newimg = image[pix]
      end

   'pclip': begin
      ndata = N_elements(image)
      newimg = sort(image)
      indx1 = long(lorej*ndata)
      indx2 = ndata - 1 - long(hirej*ndata)
      if (indx2 GE indx1) then $
       newimg = newimg[indx1:indx2]
      end

   else: begin
      newimg = image
      end
   endcase

   return, newimg
end
;------------------------------------------------------------------------------
function djs_photsky, xcen, ycen, skyrad, image, $
 salg=salg, srejalg=srejalg, smaxiter=smaxiter, $
 lorej=lorej, hirej=hirej, skyrms=skyrms, skyerr=skyerr, quick=quick, $
                      exact=exact
 
   ; Need at least 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - result = djs_photsky( xcen, ycen, skyrad, image, $'
      print, ' [ salg=, srejalg=, smaxiter=, $'
      print, ' lorej=, hirej=, skyrms=, skyerr=, /quick, /exact ] )'
      return, -1
   endif

   if (N_params() LT 2) then skyrad = 0

   if (NOT keyword_set(salg)) then begin
      if (N_params() GE 2) then salg = 'mean' $
       else salg = 'none'
   endif

   if ( (salg NE 'mean' AND salg NE 'median' AND salg NE 'mode') $
    OR NOT keyword_set(skyrad) ) then begin
      skyrms = 0.0
      skyerr = 0.0
      return, 0.0
   endif

   if ( keyword_set(srejalg) EQ 0) then srejalg = 'sigclip'
   if ( keyword_set(smaxiter) EQ 0 ) then smaxiter = 10
   case srejalg of
   'sigclip': begin
      if ( keyword_set(lorej) EQ 0 ) then lorej = 3.0
      if ( keyword_set(hirej) EQ 0 ) then hirej = 3.0
      end
   'pclip': begin
      if ( keyword_set(lorej) EQ 0 ) then lorej = 0.05
      if ( keyword_set(hirej) EQ 0 ) then hirej = 0.05
      end
   else: begin
      if ( keyword_set(lorej) EQ 0 ) then lorej = 0.0
      if ( keyword_set(hirej) EQ 0 ) then hirej = 0.0
      end
   endcase


   ; Find sky pixels and contribution from each
   xdimen = N_elements(image(*,0))
   ydimen = N_elements(image(0,*))

   if keyword_set(quick) then begin 
       quick_photfrac, xcen, ycen, skyrad, xdimen=xdimen, ydimen=ydimen, $
         pixnum=pixnum, fracs=fracs, /ragged
   endif else if keyword_set(exact) then begin 
      exact_photfrac, xcen, ycen, skyrad, xdimen=xdimen, ydimen=ydimen, $
        pixnum=pixnum, fracs=fracs
   endif else begin 
      djs_photfrac, xcen, ycen, skyrad, xdimen=xdimen, ydimen=ydimen, $
        pixnum=pixnum, fracs=fracs
   endelse 

   if (pixnum(0) NE -1) then begin

      ; Trim to only those pixels more than half filled, unless that rejects
      ; all pixels
      filled = where(fracs GT 0.5, count)
      if (count NE 0) then begin
         pixnum = pixnum[filled]
         fracs = fracs[filled]
      endif

      ; Iterate until one of the following conditions is met:
      ;   (1) the maximum number of iterations is reached
      ;   (2) no points are retained after the rejection algorithm
      ;   (3) the same sky value is returned after two iterations
      subimg = image[pixnum]
      iiter = 0
      dsky = 1.0 ; Set NE 0 to prevent stopping at first iteration
      while (iiter LE smaxiter AND N_elements(subimg) NE 0 $
       AND dsky NE 0.0) do begin

         if (iiter GT 0) then dsky = skyval

         skyval = djs_photsky_compute( subimg, salg )
         skydiff = subimg - skyval
         skyrms = sqrt( total(skydiff*skydiff) / N_elements(subimg) )
         skyerr = skyrms / sqrt( N_elements(subimg) )

         subimg = djs_photsky_reject( skyval, image[pixnum], srejalg, $
          lorej, hirej )

         if (iiter GT 0) then dsky = dsky - skyval
         iiter = iiter + 1
      endwhile

   endif else begin

      skyval = 0.0  ; No sky pixels exist in the given annulus
      skyrms = 0.0
      skyerr = 0.0

   endelse

   return, skyval
end
;------------------------------------------------------------------------------
