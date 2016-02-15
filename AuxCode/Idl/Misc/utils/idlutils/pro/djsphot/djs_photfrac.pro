;------------------------------------------------------------------------------
;+
; NAME:
;   djs_photfrac
;
; PURPOSE:
;   Create a list of pixel numbers and their fractional contribution to
;   an annular region.
;
; CALLING SEQUENCE:
;   djs_photfrac, xcen, ycen, Rvec, xdimen=, ydimen=, $
;    [ xPixNum=, yPixNum=, pixnum=, fracs=, fillfrac= ]
;
; INPUTS:
;   xcen:       X center(s)
;   ycen:       Y center(s)
;   Rvec:       Either a 2-element array with two radii to define an annulus,
;               or a scalar to define a circular aperature.
;
; OPTIONAL INPUTS:
;   xdimen:     Number of X pixels.
;   ydimen:     Number of Y pixels.
;
; OUTPUTS:
;   pixnum:     Pixel number, 0-indexed, for referencing array using one index.
;   xPixNum:    Pixel number in X, 0-indexed.
;   yPixNum:    Pixel number in Y, 0-indexed.
;   fracs:      Return value of covering fraction of the annulus
;               over the pixel number.
;   fillfrac:   Ratio of returned pixel areas to the annulus area;
;               this ratio should equal 1.0 if the aperature falls completely
;               within the image boundaries
;
; COMMENTS:
;   The total counts within this region is given by
;     totcounts = total( pData(pixnum) * fracs )
;   The area within this region is given by
;     area = total(fracs)
;   The average counts is given by
;     totcounts = total( pData(pixnum) * fracs ) / total(fracs)
;   To test for bad pixels, e.g. values greater than vmax, within
;   the aperature,
;     if (where(pData(pixnum) GT vmax) EQ -1) then <no bad values> $
;     else <bad values exist>
;
;   If no pixels within the given annulus are found, then return pixnum=-1.
;
; BUGS:
;   - can wrap around on edge of you use PixNum.  XPixNum,YPixNum do
;     not exhibit this problem
;
; PROCEDURES CALLED:
;   djs_ceil()
;   djs_floor()
;
; REVISION HISTORY:
;   Written D. Schlegel, 27 November 1996, Durham
;   Bug identified - 2 Nov 2000, D. Finkbeiner 
;-
;------------------------------------------------------------------------------
; INTERNAL SUPPORT PROCEDURES:
;
; djs_photfrac_intcirc
;------------------------------------------------------------------------------
; Return the integral of the circle with radius Radius within the pixel
; defined by x=[xA,xB] and y=[yA,yB].  Add to this the area beneath the
; pixel, i.e. the value (xB-xA)*yA.
;
; This function assumes that xB > xA >= 0 and yB > yA >= 0, and that
; the circle passes through the pixel:
;   xA^2 + yA^2 <= Radius^2 <= xB^2 + yB^2

function djs_photfrac_intcirc, xA, xB, yA, yB, Radius

   xAnorm = xA / Radius
   xBnorm = xB / Radius
   yAnorm = yA / Radius
   yBnorm = yB / Radius

   gg = 0.0 * xAnorm
   pix1 = where( yBnorm GE 1.0, count )
   if (count NE 0) then $
    gg[pix1] = xAnorm[pix1]
   pix2 = where( yBnorm LT 1.0, count )
   if (count NE 0) then $
    gg[pix2] = xAnorm[pix2] > sqrt(1.0 - yBnorm[pix2]*yBnorm[pix2])

   hh = xBnorm < sqrt(1.0 - yAnorm*yAnorm)

   result = Radius*Radius * ( (gg - xAnorm) * yBnorm + (xBnorm - hh) * yAnorm $
    + 0.5 * ( hh * sqrt(1.0 - hh*hh) + asin(hh) $
            - gg * sqrt(1.0 - gg*gg) - asin(gg) ) )
 
   return, result
end
;------------------------------------------------------------------------------
pro djs_photfrac, xcen, ycen, Rvec, xdimen=xdimen, ydimen=ydimen, $
 xPixNum=xPixNum, yPixNum=yPixNum, pixnum=pixnum, $
 fracs=fracs, fillfrac=fillfrac
 
   ; Set return values in the event of an error
   pixnum = -1L
   fracs = 0

   ; Need 2 parameters
   if N_params() LT 2 then begin
      print, 'Syntax - djs_photfrac, xcen, ycen, Rvec, $'
      print, ' xdimen=, ydimen=, xPixNum=, yPixNum=, pixnum=, $'
      print, ' fracs=, fillfrac='
      return
   endif

   ; If Rvec contains one element, then use the annulus [0,Rvec],
   ; otherwise use the annulus [Rvec[0],Rvec[1]].
   if (N_elements(Rvec) EQ 1) then begin
      Radius1 = 0.0
      Radius2 = abs(Rvec[0])
   endif else begin
      Radius1 = abs(Rvec[0])
      Radius2 = abs(Rvec[1])
   endelse
   sqRadius1 = Radius1 * Radius1
   sqRadius2 = Radius2 * Radius2

   ; Limit computations to pixels in a square centered on (xCen, yCen)
   ; which completely bounds the outer circle described by Radius2.
   iStart = 0L     > djs_floor(xcen + 0.5 - Radius2)
   iEnd   = (xdimen-1) < djs_ceil (xcen - 0.5 + Radius2)
   iLen = iEnd - iStart + 1
   if (iStart GE xdimen OR iEnd LT 0 OR iLen LT 1) then begin
;      print, 'Error - No pixels in X range' ; ???
      return
   endif

   jStart = 0L     > djs_floor(ycen + 0.5 - Radius2)
   jEnd   = (ydimen-1) < djs_ceil (ycen - 0.5 + Radius2)
   jLen = jEnd - jStart + 1
   if (jStart GE ydimen OR jEnd LT 0 OR jLen LT 1) then begin
;      print, 'Error - No pixels in Y range' ; ???
      return
   endif

   ; Compute pixel corner locations relative to (xCen,yCen)
   xA = (iStart + lindgen(iLen) - 0.5) - xcen
   xB = xA + 1
   yA = (jStart + lindgen(jLen) - 0.5) - ycen
   yB = yA + 1

   ; Split the pixel across the axis y=yCen by adding one more X pixel
   ; e.g., by adding one more element to xA and xB.
   iSplit = floor(xcen - iStart + 0.5)  ; Integer X pixel number to split
   q_splitx = iSplit GE 0 AND iSplit LT iLen
   if (q_splitx) then begin
      xA = [xA, 0.0]
      xB = [xB, xB[iSplit]]
      xB[iSplit] = 0.0
   endif

   ; Split the pixel across the axis x=xCen by adding one more Y pixel
   ; e.g., by adding one more element to yA and yB.
   jSplit = floor(ycen - jStart + 0.5)  ; Integer Y pixel number to split
   q_splity = jSplit GE 0 AND jSplit LT jLen
   if (q_splity) then begin
      yA = [yA, 0.0]
      yB = [yB, yB[jSplit]]
      yB[jSplit] = 0.0
   endif

   ; Force all relative coordinates to be non-negative and reorder
   ; values such that xB>xA and yB>yA.
   xA = abs(xA)
   xB = abs(xB)
   yA = abs(yA)
   yB = abs(yB)
   xAp = xA < xB
   xBp = xA > xB
   yAp = yA < yB
   yBp = yA > yB

   ; Compute distances to lower left corner of pixel, RadiusA,
   ; and upper right corner, RadiusB.
   nPixelx = N_elements(xAp)
   nPixely = N_elements(yAp)
   nPixels = nPixelx * nPixely
   iIndx = lindgen(nPixelx,nPixely) MOD nPixelx
   jIndx = long( lindgen(nPixelx,nPixely) / nPixelx )
   sqRadiusA = reform(xAp[iIndx] * xAp[iIndx] + yAp[jIndx] * yAp[jIndx], $
    nPixelx, nPixely)
   sqRadiusB = reform(xBp[iIndx] * xBp[iIndx] + yBp[jIndx] * yBp[jIndx], $
    nPixelx, nPixely)

   ; Integrate within the annulus defined by the circles [Radius1,Radius2]
   ; within each pixel.
   Integral = fltarr(nPixelx,nPixely)

   qpix0 = ( sqRadiusB GT sqRadius1 AND sqRadiusA LT sqRadius2 )
   pix1 = where( sqRadiusB LE sqRadius2 AND qpix0, count)
   if (count NE 0) then $
    Integral[pix1] = ( xBp[iIndx[pix1]] - xAp[iIndx[pix1]] ) * yBp[jIndx[pix1]]
   pix2 = where( sqRadiusB GT sqRadius2 AND qpix0, count)
   if (count NE 0) then $
    Integral[pix2] = djs_photfrac_intcirc( xAp[iIndx[pix2]], xBp[iIndx[pix2]], $
     yAp[jIndx[pix2]], yBp[jIndx[pix2]], Radius2 )

   pix1 = where( sqRadiusA GE sqRadius1 AND qpix0, count)
   if (count NE 0) then $
    Integral[pix1] = Integral[pix1] $
     - ( xBp[iIndx[pix1]] - xAp[iIndx[pix1]] ) * yAp[jIndx[pix1]]
   pix2 = where( sqRadiusA LT sqRadius1 AND qpix0, count)
   if (count NE 0) then $
    Integral[pix2] = Integral[pix2] $
     - djs_photfrac_intcirc( xAp[iIndx[pix2]], xBp[iIndx[pix2]], $
     yAp[jIndx[pix2]], yBp[jIndx[pix2]], Radius1 )

   ; Collapse the split pixels back into the original pixels
   if (q_splity) then $
    Integral[0:nPixelx-1,jSplit] = Integral[0:nPixelx-1,jSplit] $
     + Integral[0:nPixelx-1,nPixely-1]
   if (q_splitx) then $
    Integral[iSplit,0:nPixely-1] = Integral[iSplit,0:nPixely-1] $
     + Integral[nPixelx-1,0:nPixely-1]

   ; Set the return values
   xPixNum = iIndx[0:nPixelx-1-q_splitx,0:nPixely-1-q_splity] + iStart
   yPixNum = jIndx[0:nPixelx-1-q_splitx,0:nPixely-1-q_splity] + jStart
   fracs = Integral[0:nPixelx-1-q_splitx,0:nPixely-1-q_splity]

   ; Limit the return values to only those pixels with non-zero contributions
   pix1 = where(fracs NE 0.0)
   if (pix1[0] EQ -1) then return ; ??? Should never happen!
   xPixNum = xPixNum[pix1]
   yPixNum = yPixNum[pix1]
   pixnum = 0L + xPixNum + xdimen * yPixNum
   fracs = fracs[pix1]

   ; Test to see if aperature exceeds image boundary by computing the 
   ; ratio of the filled pixels to the area of the annulus.  If all
   ; pixels are within the image boundary, then fillfrac=1.0.
   fillfrac = total(fracs) / (!pi * (sqRadius2 - sqRadius1))

   return
end
;------------------------------------------------------------------------------
