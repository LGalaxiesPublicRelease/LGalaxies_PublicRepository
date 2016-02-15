;------------------------------------------------------------------------------
;+
; NAME:
;   djs_photcen
;
; PURPOSE:
;   Recenter an object position.
;
; CALLING SEQUENCE:
;   djs_photcen, xcen, ycen, image, $
;    [ calg=, cbox=, cmaxiter=, cmaxshift=, fwhm=, /fixfw, ceps=, qmaxshift= ]
;
; INPUTS:
;   xcen:       X center(s)
;   ycen:       Y center(s)
;   image:      2-dimensional image
;
; OPTIONAL INPUTS:
;   calg:       Centering algorithm.  Choose from iweight, gauss1, gauss2, none.
;                 iweight = intensity-weighted center, computed independently
;                           in both X and Y
;                 gauss1  = gaussian fit, including a constant term, computed
;                           independently in both X and Y
;                 gauss2  = 2D gaussian fit, including a constant term
;                 none    = no centering
;               Default to iweight.
;   cbox:       Centering box width.  Default to 7.
;   cmaxiter:   Maximum number of iterations for centering algorithm.
;               Default to 10.
;   cmaxshift:  Maximum center shift.  If this shift is exceeded in either
;               X or Y, then return the center position passed in XCEN, YCEN.
;               A value of 0 imposes no maximum shift.  Default to 3.
;   fwhm:       FWHM for gauss1 or gauss2 centering algorithms.  If this is
;               a scalar, then the same value is used for both X and Y.
;               If this is a vector, then the first two elements are used
;               for X and Y respectively.
;   fixfw:      If set and nonzero, then fix the FWHM for gauss1 or gauss2 fits.
;   ceps:       Stop iterating when relative shifts in both X and Y are less
;               than this value.  Default to 0.
;
; OUTPUTS:
;   xcen:       Re-centered X position (modified).
;   ycen:       Re-centered X position (modified).
;
; OPTIONAL OUTPUS:
;   qmaxshift:  Return 1 if maximum shift was reached in either X or Y.
;               Return 0 otherwise.
;
; PROCEDURES CALLED:
;   curvefit()
;   djs_ceil()
;   djs_floor()
;   gauss2dfit()
;
; REVISION HISTORY:
;   01-Dec-1996  Written by D. Schlegel, Durham.
;   10-Aug-1998  Added option for calg='gauss2' (DJS).
;   01-Jun-2000  Major revisions: change XYCEN to XCEN,YCEN; remove use
;                of FITS headers; make IDL 5 complient (DJS).
;-
;------------------------------------------------------------------------------
; INTERNAL SUPPORT PROCEDURES:
;
; djs_photcen_gfunc1
; djs_photcen_gfit1()
;------------------------------------------------------------------------------
; Function of the form
;   F(x) = a * exp(-e^2) + c
; where e=(x-x0)/sigma, and the fitting parameters avec=[a,x0,c,sigma].
; This function can be passed to curvefit().

pro djs_photcen_gfunc1, xvec, avec, Fval, pderiv
   common djs_photcen_com, fwhcom

   if (N_elements(avec) LE 3) then fixfw = 1 $
    else fixfw = 0

   if (fixfw EQ 0) then fwhcom = avec[3]

   ee = (xvec - avec[1]) / fwhcom
   ee2 = ee * ee
   bx = exp(-0.5*ee2)
   Fval = avec[0] * bx + avec[2]

   if (N_params() GE 4) then begin
      ; Calculate partial derivatives
      dFda = bx
      dFdx0 = bx * ee / fwhcom
      dFdc = 0*bx + 1.0
      pderiv = [[bx], [dFdx0], [dFdc]]
      if (fixfw EQ 0) then begin ;  FWHM is not fixed, so include it in fit
         dFdsig = dFdx0 * ee
         pderiv = [[pderiv], [dFdsig]]
      endif
   endif

end

;------------------------------------------------------------------------------
function djs_photcen_gfit1, xvec, yvec, wvec, fwhm, xcen, fixfw

   common djs_photcen_com, fwhcom

   ; If too few data points, then return with no change to xcen
   if (N_elements(xvec) LT 4) then begin
      print, 'Error - Too few data points for gaussian fit'
      return, xcen
   endif

   ; Determine initial guesses
   c0 = min(yvec)
   sig0 = fwhm / 2.354820
   a0 = max(yvec, imax) - c0
   x0 = xvec[imax]
   avec = [a0, x0, c0]

   if (fixfw NE 0) then fwhcom = fwhm $
    else avec = [avec, sig0]

   yfit = curvefit( xvec, yvec, wvec, avec, $
    function_name='djs_photcen_gfunc1' )

   newx = avec[1]
   return, newx
end

;------------------------------------------------------------------------------
pro djs_photcen, xcen, ycen, image, $
 calg=calg, cbox=cbox, cmaxiter=cmaxiter, cmaxshift=cmaxshift, $
 fwhm=fwhm, fixfw=fixfw, ceps=ceps, qmaxshift=qmaxshift

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - result = djs_photcen( xcen, ycen, image, $'
      print, ' [ calg=, cbox=, cmaxiter=, cmaxshift=, $'
      print, ' fwhm=, /fixfw, ceps= ] )'
      return
   endif

   ;----------
   ; If XCEN and YCEN are arrays, then call this routine recursively

   nx = n_elements(xcen)
   ny = n_elements(ycen)
   if (nx NE ny) then $
    message, 'Dimensions of NX and NY do not agree'
;   if (nx GT 1) then begin
   if (size(xcen))[0] NE 0 then begin
      qmaxshift = bytarr(nx)
      for i=0, nx-1 do begin
         xtmp = xcen[i]
         ytmp = ycen[i]
         djs_photcen, xtmp, ytmp, image, $
          calg=calg, cbox=cbox, cmaxiter=cmaxiter, cmaxshift=cmaxshift, $
          fwhm=fwhm, fixfw=fixfw, ceps=ceps, qmaxshift=qtmp
         xcen[i] = xtmp
         ycen[i] = ytmp
         qmaxshift[i] = qtmp
      endfor
      return
   endif

   if ( keyword_set(calg) EQ 0 ) then calg = 'iweight'
   if ( keyword_set(cbox) EQ 0 ) then cbox = 7
   if ( keyword_set(cmaxiter) EQ 0 ) then cmaxiter = 10
   if ( keyword_set(cmaxshift) EQ 0 ) then cmaxshift = 3.0
   if ( keyword_set(ceps) EQ 0 ) then ceps = 0.0

   if ( keyword_set(fwhm) EQ 0 ) then $
      fwhvec = [1.0, 1.0] $
   else if ( keyword_set(fwhm) EQ 1 ) then $
      fwhvec = [fwhm, fwhm] $
   else $
      fwhvec = [fwhm(0), fwhm(1)]
   if ( keyword_set(fixfw) EQ 0 ) then fixfw = 0

   ; Return if no centering is to be performed
   if (calg EQ 'none') then begin
      return
   endif

   ; Use the data array dimensions
   dims = size(image, /dimens)
   naxis1 = dims[0]
   naxis2 = dims[1]

   Radius = 0.5 * cbox

   ; Iterate until one of the following conditions is met:
   ;   (1) the maximum number of iterations is reached
   ;   (2) no pixels are within the centering box
   ;   (3) the maximum shift in both X and Y has been exceeded
   ;   (4) the same center position is returned after two iterations,
   ;       to within ceps of each other in both X and Y
   iiter = 0
   dcen = 2*ceps + 1.0 ; Set > ceps to prevent stopping at first iteration
   qmaxshift = 0B
   while (iiter LE cmaxiter AND qmaxshift EQ 0 $
    AND max(abs(dcen)) GT ceps) do begin

      if (iiter GT 0) then dcen = [xcen, ycen]

      ; Limit computations to pixels in a square centered on (xCen, yCen)
      ; which bounds the box with radius Radius.

      ; Trim the X radius to be the largest symmetric box that falls within
      ; the image boundaries.
      xRad = min([Radius, (xcen > 0), ((naxis1-xcen) > 0)])
      iStart = djs_floor(xcen + 0.5 - xRad)
      iEnd   = djs_ceil(xcen - 0.5 + xRad)
      if (iStart GE naxis1 OR iEnd LT 0) then begin
         print, 'Error - No pixels in X range'
         return
      endif
      iLen = iEnd - iStart + 1
 
      ; Trim the Y radius to be the largest symmetric box that falls within
      ; the image boundaries.
      yRad = min([Radius, (ycen > 0), ((naxis2-ycen) > 0)])
      jStart = djs_floor(ycen + 0.5 - yRad)
      jEnd   = djs_ceil(ycen - 0.5 + yRad)
      if (jStart GE naxis2 OR jEnd LT 0) then begin
         print, 'Error - No pixels in Y range'
         return
      endif
      jLen = jEnd - jStart + 1

      ; Compute pixel corner locations relative to (xCen,yCen)
      xA = iStart + indgen(iLen) - 0.5 - xcen
      xB = xA + 1
      yA = jStart + indgen(jLen) - 0.5 - ycen
      yB = yA + 1

      ; Trim pixel sizes to only those coordinates falling within the
      ; centering box
      xA = xA > (-xRad)
      xB = xB < xRad
      yA = yA > (-yRad)
      yB = yB < yRad

      ; Determine the indices for the pixel
      nPixelx = N_elements(xA)
      nPixely = N_elements(yA)
      nPixels = nPixelx * nPixely
      iIndx = lindgen(nPixelx,nPixely) MOD nPixelx
      jIndx = long( lindgen(nPixelx,nPixely) / nPixelx ) + 0L
      xPixNum = iIndx + iStart
      yPixNum = jIndx + jStart
      pixnum = 0L + xPixNum + naxis1 * yPixNum

      ; Compute contribution of each pixel
      fracs = (xB[iIndx[*]] - xA[iIndx[*]]) * (yB[jIndx[*]] - yA[jIndx[*]])

      ; Create a temporary data array with the box of pixels
      subimg = image[pixnum]

      case calg of
      'iweight': begin
         ; Intensity-weighted centering
         norm = total(subimg[*] * fracs)

         ; Insist that the total flux is positive
         if (norm GT 0) then begin
            ; Work with xPixNum-xcen and yPixNum-ycen for numerical stability
            newi = total( subimg[*] * (xPixNum[*]-xcen) * fracs) / norm + xcen
            newj = total( subimg[*] * (yPixNum[*]-ycen) * fracs) / norm + ycen
            xcen = newi
            ycen = newj
         endif

         end

      'gauss1': begin
         ; One-dimensional gaussian fit, indepently fit in X and Y
         ; Collapse the sub-image into two one-dimensional arrays

         ; Fit X center
         pTemp = total(subimg, 2)
         fracx = xB - xA
         newi = djs_photcen_gfit1( xPixNum[*,0], pTemp, fracx, $
          fwhvec[0], xcen, fixfw )
         xcen = newi

         ; Fit Y center
         pTemp = total(subimg, 1)
         fracy = yB - yA
         newj = djs_photcen_gfit1( transpose(yPixNum[0,*]), pTemp, fracy, $
          fwhvec[1], ycen, fixfw )
         ycen = newj

         end

      'gauss2': begin
         ; Two-dimensional gaussian fit, simultaneously fitting in X and Y

         gres = gauss2dfit( subimg, gcoeff, /tilt )
         xcen = gcoeff[4]+iStart
         ycen = gcoeff[5]+jStart

         end

      'none': begin
         end

      else: begin
         print, 'Error - Invalid centering algorithm'
         return
         end
      endcase

      ; Test to see if maximum shift was reached in either X or Y.
      ; If so, then reset coordinates to original center and set the
      ; qmaxshift flag.
      if (keyword_set(cmaxshift)) then begin
         xold = xcen
         yold = ycen
         xcen = (xcen < xcen + cmaxshift) > xcen - cmaxshift
         ycen = (ycen < ycen + cmaxshift) > ycen - cmaxshift
         if (xold NE xcen OR yold NE ycen) then qmaxshift = 1B
      endif

      if (iiter GT 0) then dcen = dcen - [xcen,ycen]
      iiter = iiter + 1

   endwhile

   return
end
;------------------------------------------------------------------------------
