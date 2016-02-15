pro xy2ad, x, y, astr, a, d
;+
; NAME:
;     XY2AD
;
; PURPOSE:
;     Compute R.A. and Dec from X and Y and a FITS astrometry structure
; EXPLANATION:
;     The astrometry structure must first be extracted by EXTAST from a FITS
;     header.   The offset from the reference pixel is computed and the CD 
;     matrix is applied.     Of distortion is present then this is corrected.
;     If a WCS projection (Calabretta & Greisen 2002, A&A, 395, 1077) is 
;     present, then the procedure WCSXY2SPH is used to compute astronomical
;     coordinates.    Angles are returned in  degrees.
;   
;     XY2AD is meant to be used internal to other procedures.  
;     For interactive purposes use XYAD.
;
; CALLING SEQUENCE:
;     XY2AD, x, y, astr, a, d   
;
; INPUTS:
;     X     - row position in pixels, scalar or vector
;     Y     - column position in pixels, scalar or vector
;           X and Y should be in the standard IDL convention (first pixel is
;           0), and not the FITS convention (first pixel is 1). 
;     ASTR - astrometry structure, output from EXTAST procedure containing:
;        .CD   -  2 x 2 array containing the astrometry parameters CD1_1 CD1_2
;               in DEGREES/PIXEL                                   CD2_1 CD2_2
;        .CDELT - 2 element vector giving physical increment at reference pixel
;        .CRPIX - 2 element vector giving X and Y coordinates of reference pixel
;               (def = NAXIS/2)
;        .CRVAL - 2 element vector giving R.A. and DEC of reference pixel 
;               in DEGREES
;        .CTYPE - 2 element vector giving projection types 
;        .LONGPOLE - scalar longitude of north pole
;        .LATPOLE - scalar giving native latitude of the celestial pole
;        .PV2 - Vector of projection parameter associated with latitude axis
;             PV2 will have up to 21 elements for the ZPN projection, up to 3
;             for the SIN projection and no more than 2 for any other
;             projection
;        .DISTORT - Optional substructure specifying distortion parameters
;                  
;
; OUTPUT:
;     A - R.A. in DEGREES, same number of elements as X and Y
;     D - Dec. in DEGREES, same number of elements as X and Y
;
; RESTRICTIONS:
;       Note that all angles are in degrees, including CD and CRVAL
;       Also note that the CRPIX keyword assumes an FORTRAN type
;       array beginning at (1,1), while X and Y give the IDL position
;       beginning at (0,0).   No parameter checking is performed.
;
; NOTES:
;      AD2XY tests for presence of WCS coordinates by the presence of a dash 
;      in the 5th character position in the value of CTYPE (e.g 'DEC--SIN').       
; PROCEDURES USED:
;       TAG_EXIST(), WCSXY2SPH
; REVISION HISTORY:
;       Written by R. Cornett, SASC Tech., 4/7/86
;       Converted to IDL by B. Boothman, SASC Tech., 4/21/86
;       Perform CD  multiplication in degrees  W. Landsman   Dec 1994
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Understand reversed X,Y (X-Dec, Y-RA) axes,   W. Landsman  October 1998
;       Consistent conversion between CROTA and CD matrix W. Landsman Oct. 2000
;       No special case for tangent projection W. Landsman June 2003
;       Work for non-WCS coordinate transformations W. Landsman Oct 2004
;- 
 compile_opt idl2
 if N_params() LT 4 then begin
        print,'Syntax -- XY2AD, x, y, astr, a, d'
        return
 endif
 radeg = 180.0d/!DPI                  ;Double precision !RADEG

 cd = astr.cd
 crpix = astr.crpix
 cdelt = astr.cdelt
 if cdelt[0] NE 1.0 then begin 
         cd[0,0] = cd[0,0]*cdelt[0] & cd[0,1] = cd[0,1]*cdelt[0]
         cd[1,1] = cd[1,1]*cdelt[1] & cd[1,0] = cd[1,0]*cdelt[1]
  endif


 xdif = x - (crpix[0]-1)            
 ydif = y - (crpix[1]-1)
 
 if tag_exist(astr,'DISTORT') then begin
      if astr.distort.name EQ 'SIP' then begin
           distort  = astr.distort
           a = distort.a
           b = distort.b
           na = ((size(a,/dimen))[0])
           xdif1 = xdif
           ydif1 = ydif
           
           for i=0,na-1 do begin
               for j=0,na-1 do begin
                  if a[i,j] NE 0.0 then xdif1 = xdif1 + xdif^i*ydif^j*a[i,j]            
                  if b[i,j] NE 0.0 then ydif1 = ydif1 + xdif^i*ydif^j*b[i,j]
           endfor
           endfor

           xdif = xdif1
           ydif = ydif1
           
      endif
 endif

 xsi = cd[0,0]*xdif + cd[0,1]*ydif   ;Can't use matrix notation, in
 eta = cd[1,0]*xdif + cd[1,1]*ydif   ;case X and Y are vectors

 ctype = astr.ctype
 crval = astr.crval
 coord = strmid(ctype,0,4)
 reverse = ((coord[0] EQ 'DEC-') and (coord[1] EQ 'RA--')) or $
           ((coord[0] EQ 'GLAT') and (coord[1] EQ 'GLON')) or $
           ((coord[0] EQ 'ELAT') and (coord[1] EQ 'ELON'))

 if reverse then begin
     crval = rotate(crval,2)
     temp = xsi & xsi = eta & eta = temp
 endif

 if strmid(ctype[0],4,1) EQ '-' then begin
 WCSXY2SPH, xsi, eta, a, d, CTYPE = ctype[0:1], PV2 = astr.pv2, $
        LONGPOLE = astr.longpole, CRVAL = crval, LATPOLE = astr.latpole
 endif else begin
         a = xsi & d = eta
 endelse
 return
 end
