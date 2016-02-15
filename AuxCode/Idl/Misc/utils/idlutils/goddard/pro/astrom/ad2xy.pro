pro ad2xy, a, d, astr, x, y
;+
; NAME:
;     AD2XY
; PURPOSE:
;     Compute X and Y from native coordinates and a FITS  astrometry structure
; EXPLANATION:
;     If a WCS projection (Calabretta & Greisen 2002, A&A, 395, 1077) is 
;     present, then the procedure WCSXY2SPH is used to compute native 
;     coordinates.   If distortion is present then this is corrected.  
;     In all cases, the inverse of the CD matrix is applied and offset 
;     from the reference pixel to obtain X and Y. 
;
;     AD2XY is generally meant to be used internal to other procedures.   For 
;     interactive purposes, use ADXY.
;
; CALLING SEQUENCE:
;     AD2XY, a ,d, astr, x, y   
;
; INPUTS:
;     A -     R.A. or longitude in DEGREES, scalar or vector
;     D -     Dec. or longitude in DEGREES, scalar or vector
;     ASTR - astrometry structure, output from EXTAST procedure containing:
;        .CD   -  2 x 2 array containing the astrometry parameters CD1_1 CD1_2
;               in DEGREES/PIXEL                                   CD2_1 CD2_2
;        .CDELT - 2 element vector giving increment at reference point in
;               DEGREES/PIXEL
;        .CRPIX - 2 element vector giving X and Y coordinates of reference pixel
;               (def = NAXIS/2) in FITS convention (first pixel is 1,1)
;        .CRVAL - 2 element vector giving coordinates of the reference pixel 
;               in DEGREES
;        .CTYPE - 2 element vector giving projection types 
;        .LONGPOLE - scalar longitude of north pole (default = 180) 
;        .PV2 - Vector of additional parameter (e.g. PV2_1, PV2_2) needed in 
;               some projections
;        .DISTORT - Optional substructure specifying distortion parameters
;
; OUTPUTS:
;     X     - row position in pixels, scalar or vector
;     Y     - column position in pixels, scalar or vector
;
;     X,Y will be in the standard IDL convention (first pixel is 0), and
;     *not* the FITS convention (first pixel is 1)
; NOTES:
;      AD2XY tests for presence of WCS coordinates by the presence of a dash 
;      in the 5th character position in the value of CTYPE (e.g 'DEC--SIN').       
; PROCEDURES USED:
;       TAG_EXIST(), WCSSPH2XY
; REVISION HISTORY:
;     Converted to IDL by B. Boothman, SASC Tech, 4/21/86
;     Use astrometry structure,  W. Landsman      Jan. 1994   
;     Do computation correctly in degrees  W. Landsman       Dec. 1994
;     Only pass 2 CRVAL values to WCSSPH2XY   W. Landsman      June 1995
;     Don't subscript CTYPE      W. Landsman       August 1995        
;     Converted to IDL V5.0   W. Landsman   September 1997
;     Understand reversed X,Y (X-Dec, Y-RA) axes,   W. Landsman  October 1998
;     Consistent conversion between CROTA and CD matrix W. Landsman October 2000
;     No special case for tangent projection W. Landsman June 2003
;     Work for non-WCS coordinate transformations W. Landsman Oct 2004
;-
 On_error,2
 compile_opt idl2

 if N_params() lT 4 then begin
        print,'Syntax -- AD2XY, a, d, astr, x, y'
        return
 endif

 radeg = 180.0D/!DPI                 ;Double precision !RADEG
 ctype = astr.ctype
 crval = astr.crval

 coord = strmid(ctype,0,4)
 reverse = ((coord[0] EQ 'DEC-') and (coord[1] EQ 'RA--')) or $
           ((coord[0] EQ 'GLAT') and (coord[1] EQ 'GLON')) or $
           ((coord[0] EQ 'ELAT') and (coord[1] EQ 'ELON'))
 if reverse then crval = rotate(crval,2)        ;Invert CRVAL?

 if (ctype[0] EQ '') then begin   
      ctype = ['RA---TAN','DEC--TAN']
      message,'No CTYPE specified - assuming TANgent projection',/INF
 endif      
     
  spherical = strmid(astr.ctype[0],4,1) EQ '-'
  if spherical then begin
  wcssph2xy, a, d, xsi, eta, CTYPE = ctype, PV2 = astr.pv2, $
        LONGPOLE = astr.longpole, CRVAL = crval, LATPOLE = astr.latpole
  endif else begin
        xsi = a & eta = d
  endelse	
  cd = astr.cd
  cdelt = astr.cdelt

  if cdelt[0] NE 1.0 then begin
         cd[0,0] = cd[0,0]*cdelt[0] & cd[0,1] = cd[0,1]*cdelt[0]
         cd[1,1] = cd[1,1]*cdelt[1] & cd[1,0] = cd[1,0]*cdelt[1]
     endif

 if reverse then begin
     temp = xsi &  xsi = eta & eta = temp
 endif

 crpix = astr.crpix - 1
 cdinv = invert(cd)
 xdif = ( cdinv[0,0]*xsi + cdinv[0,1]*eta  )
 ydif = ( cdinv[1,0]*xsi + cdinv[1,1]*eta  )

 if tag_exist(astr,'DISTORT') then begin
      if astr.distort.name EQ 'SIP' then begin
           distort  = astr.distort
           ap = distort.ap
           bp = distort.bp
           na = ((size(ap,/dimen))[0])
           xdif1 = xdif
           ydif1 = ydif
           
           for i=0,na-1 do begin
               for j=0,na-1 do begin
                  if ap[i,j] NE 0.0 then xdif1 = xdif1 + xdif^i*ydif^j*ap[i,j]            
                  if bp[i,j] NE 0.0 then ydif1 = ydif1 + xdif^i*ydif^j*bp[i,j]
           endfor
           endfor

           xdif = xdif1
           ydif = ydif1
           
      endif
 endif

 x = xdif + crpix[0] 
 y = ydif + crpix[1] 
 return
 end
