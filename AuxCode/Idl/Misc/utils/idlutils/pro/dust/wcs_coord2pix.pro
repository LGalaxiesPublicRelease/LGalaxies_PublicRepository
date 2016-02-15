;------------------------------------------------------------------------------
;+
; NAME:
;   wcs_coord2pix
;
; PURPOSE:
;   Transform from a world coordinate system (WCS) to (x,y) pixel numbers.
;   This function returns the ZERO-INDEXED pixel position.
;
;   If the FRACTIONAL flag is set, then a fractional pixel position is
;   returned.
;
; CALLING SEQUENCE:
;   wcs_coord2pix, lonvec, latvec, hdr, xpix, ypix, [ /fractional ]
;
; INPUTS:
;   lonvec:     Galactic longitude [degrees]
;   latvec:     Galactic latitude [degrees]
;   hdr:        FITS header data
;
; OPTIONAL INPUTS:
;   fractional: Set this to return fractional rather than integer
;               pixel numbers
;
; OUTPUTS:
;   xpix:       X pixel position in range [0,naxis1-1]
;   ypix:       Y pixel position in range [0,naxis2-1]
;
; COMMENTS:
;   This conversion routine is based upon the FITS standard of Griesen &
;   Calabretta (1996).  Equation numbers refer to that paper.
;
;   At present, only zenithal equal area (ZEA) projections are supported.
;
;   A special case is included for the Lambert projection as defined by
;   Schlegel, Finkbeiner & Davis (1998).  Galactic latitude runs clockwise
;   from X-axis for NGP, counterclockwise for SGP.
;   (CRPIX1, CRPIX2) define the 1-indexed central pixel location of
;   the pole.  For example, if CRPIX1=512, then the pole is exactly in
;   the middle of 1-indexed pixel number 512; if CRPIX1=512.5, then
;   the pole falls between pixel numbers 512 and 513.
;
; PROCEDURES CALLED:
;   djs_angpos()
;   sxpar()
;
; REVISION HISTORY:
;   Written by D. Schlegel, 30 May 1996, Durham
;   18-Mar-1999  Use of TEMPORARY function to improve memory use
;                (D. Finkbeiner).
;   30-Mar-1999  Re-written as a general forward conversion for world
;                coordinate systems (WCS), though most not yet implemented
;                (DJS, DPF).
;-
;------------------------------------------------------------------------------
pro wcs_coord2pix, lonvec, latvec, hdr, xpix, ypix, fractional=fractional

   ; Need five parameters
   if (N_params() LT 5) then begin
      print, 'Syntax - wcs_coord2pix, lonvec, latvec, hdr, xpix, ypix, /fractional'
      return
   endif

   ; Test hdr
   if (keyword_set(hdr) EQ 0) then begin 
      print, 'wcs_coord2pix - ERROR:  hdr must contain FITS header'
      return
   endif 

   szhead = size(hdr)
   if (szhead[0] NE 1) OR (szhead[2] NE 7) then begin 
      print, 'wcs_coord2pix - ERROR:  hdr must contain FITS header'
      return
   endif 

   ; Determine the projection
   ctype1 = sxpar(hdr, 'CTYPE1')
   ctype2 = sxpar(hdr, 'CTYPE2')
   ctype1_coord = (str_sep(ctype1,'-'))[0]
   ctype2_coord = (str_sep(ctype2,'-'))[0]
   ctype1_proj = strmid(ctype1,5, 255) ; for compatibility with older IDL
   ctype2_proj = strmid(ctype2,5, 255)

   ; Get information from FITS header necessary for any projection
   crval1 = sxpar(hdr, 'CRVAL1')
   crval2 = sxpar(hdr, 'CRVAL2')
   crpix1 = sxpar(hdr, 'CRPIX1')
   crpix2 = sxpar(hdr, 'CRPIX2')
   lonpole = sxpar(hdr, 'LONPOLE', count=ct1)
   if (ct1 EQ 0) then lonpole = sxpar(hdr, 'LONGPOLE')

   cd1_1 = sxpar(hdr, 'CDELT1', count=ct1)
   cd2_2 = sxpar(hdr, 'CDELT2', count=ct2)
   if (ct1 EQ 1 AND ct2 EQ 1) then begin
      cd1_2 = 0.0
      cd2_1 = 0.0
   endif else begin
      cd1_1 = sxpar(hdr, 'CD1_1', count=ct1)
      cd1_2 = sxpar(hdr, 'CD1_2', count=ct2)
      cd2_1 = sxpar(hdr, 'CD2_1', count=ct3)
      cd2_2 = sxpar(hdr, 'CD2_2', count=ct4)
      if (ct1+ct2+ct3+ct4 EQ 0) then $
       message, 'No CDELT cards or CD matrix in header'
   endelse

   ; Set the conversion from degrees to radians to the same precision
   ; as the input coordinates
   if (size(lonvec,/tname) EQ 'DOUBLE' OR size(latvec,/tname) EQ 'DOUBLE') $
    then begin
      dradeg = 180.d0 / !dpi
      sqtwo = sqrt(2.d0)
   endif else begin
      dradeg = 180. / !pi
      sqtwo = sqrt(2.)
   endelse

   if (ctype1 EQ 'LAMBERT--X' AND ctype2 EQ 'LAMBERT--Y') then begin

      nsgp = sxpar(hdr, 'LAM_NSGP')
      scale = sxpar(hdr, 'LAM_SCAL')

;      rho = float( sqrt( 1. - nsgp * sin(double(latvec)/dradeg) ) )
      rho = sqtwo * sin((45. - 0.5 * nsgp * latvec)/dradeg)
      xpix =         scale * rho * cos(lonvec/dradeg)
      ypix = -nsgp * scale * rho * sin(lonvec/dradeg)

      xpix = temporary(xpix) + (crpix1 - crval1 - 1.0)
      ypix = temporary(ypix) + (crpix2 - crval2 - 1.0)

   endif else if (ctype1_proj EQ 'ZEA') then begin

      if (ctype1_proj NE ctype2_proj) then $
       message, 'Inconsistent projection names in CTYPE keywords'

      ; ROTATION
      ; Equn (4) - degenerate cases only
      if (crval2 EQ 90) then begin
         theta = latvec
         phi = djs_angpos( lonvec + (180 + lonpole - crval1) )
      endif else if (crval2 EQ -90) then begin
         theta = -latvec
         phi = djs_angpos( (lonpole + crval1) - lonvec)
      endif else begin
         message, 'This value of CRVAL2 not yet supported.'
      endelse

      ; FORWARD MAP PROJECTION
      ; Equn (26)
      Rtheta = 2 * dradeg * sin((0.5 / dradeg) * (90 - temporary(theta)))

      ; Equns (10), (11)
      xtemp = Rtheta * sin(phi / dradeg)
      ytemp = - temporary(Rtheta) * cos(temporary(phi) / dradeg)

      ; SCALE FROM PHYSICAL UNITS
      ; Equn (3) after inverting the matrix
      if (cd1_1 EQ 0 AND cd2_2 EQ 0) then begin
         xpix = temporary(xtemp) / cd1_1 + (crpix1 - 1.0)
         ypix = temporary(ytemp) / cd2_2 + (crpix2 - 1.0)
      endif else begin
         denom = cd1_1 * cd2_2 - cd1_2 * cd2_1
         xpix = (cd2_2 * xtemp - cd1_2 * ytemp) / denom + (crpix1 - 1.0)
         ypix = (cd1_1 * ytemp - cd2_1 * xtemp) / denom + (crpix2 - 1.0)
xtemp = 0
ytemp = 0
      endelse

   endif else begin

      message, 'Unsupported projection method in CTYPE keywords'

   endelse

   ; Round to nearest pixel if /FRACTIONAL keyword not set.
   ; Note that this may return pixels off the edge, depending on
   ; numerical details. 
   if (NOT keyword_set(fractional)) then begin
      xpix = fix(temporary(xpix) + 0.5)
      ypix = fix(temporary(ypix) + 0.5)
   endif

   return
end
;------------------------------------------------------------------------------
