pro gsssputast, hdr, astr
;+
; NAME:
;    GSSSPUTAST
; PURPOSE:
;    Put GSSS astrometry parameters into a given FITS header.
;
; CALLING SEQUENCE:
;     gsssputast, hdr, astr
;
; INPUTS:
;     hdr   - FITS header, string array.   HDR will be updated to contain
;               the supplied astrometry.
;     astr  - IDL structure containing the GSSS Astrometry information
;               .CTYPE[2]   =  ['RA---GSS','DEC--GSS'] 
;               .CRVAL[2]   = plate center Ra, Dec (from PLTRAH, PLTRAM etc.)
;               .XLL,.YLL   = offsets lower lefthand corner (integers)
;               .AMDX[13],.AMDY[13] = 12 transformation coefficients
;               .XSZ,.YSZ   = X and Y pixel size in microns
;               .PLTSCL     = plate scale in arc sec/mm
;               .PPO3,.PPO6 = orientation coefficients
;
; OUTPUTS:
;      hdr - (Modified) FITS header with updated astrometry cards.
;            A brief HISTORY record is also added.
;
; REVISION HISTORY:
;   Written 4 Nov 2000 - D. Finkbeiner
;-
;------------------------------------------------------------------------------
  npar = N_params()
   
  if (npar LT 2) then begin ; Was header supplied?
     print, 'Syntax: GSSSPUTAST, hdr, astr '
     return
  end

  IF NOT (size(astr,/tname) EQ 'STRUCT') THEN BEGIN 
     print, 'Syntax: GSSSPUTAST, hdr, astr '
     return
  end

  ctype = astr.ctype

  IF (ctype[0] NE 'RA---GSS') OR (ctype[1] NE 'DEC--GSS') THEN BEGIN 
     print, 'GSSSPUTAST: Header must be type GSS'
     return
  ENDIF 

  sxaddpar, hdr, 'CTYPE1',   ctype[0], ' Coordinate Type'
  sxaddpar, hdr, 'CTYPE2',   ctype[1], ' Coordinate Type'
  sxaddpar, hdr, 'CRPIX1',   astr.xll, ' Ref Pixel X'
  sxaddpar, hdr, 'CRPIX2',   astr.yll, ' Ref Pixel Y'
  sxaddpar, hdr, 'XPIXELSZ', astr.xsz, ' Pixel size'
  sxaddpar, hdr, 'YPIXELSZ', astr.ysz, ' Pixel size'
  sxaddpar, hdr, 'PPO1',     0.,       ' Placeholder, so extast works'
  sxaddpar, hdr, 'PPO3',     astr.ppo3, ''
  sxaddpar, hdr, 'PPO6',     astr.ppo6, ''

  sxaddpar, hdr, 'PLTSCALE', astr.pltscl, ' Plate scale'

  RA = astr.crval[0]/15.
  pltrah = fix(RA)
  pltram = fix((RA-pltrah)*60.)
  pltras = ((RA-pltrah)*60.-pltram)*60.

  dec = abs(astr.crval[1])
  pltdecd = fix(dec)
  pltdecm = fix((dec-pltdecd)*60.)
  pltdecs = ((dec-pltdecd)*60.-pltdecm)*60.
  pltdecsn = (astr.crval[1] GT 0.) ? '+' : '-'

  sxaddpar, hdr, 'PLTRAH',   pltrah   
  sxaddpar, hdr, 'PLTRAM',   pltram  
  sxaddpar, hdr, 'PLTRAS',   pltras  
  sxaddpar, hdr, 'PLTDECSN', pltdecsn
  sxaddpar, hdr, 'PLTDECD',  pltdecd 
  sxaddpar, hdr, 'PLTDECM',  pltdecm 
  sxaddpar, hdr, 'PLTDECS',  pltdecs 

; Now do the 13 amdx coeffs

  amdx_str = strarr(13)
  amdy_str = strarr(13)

  amdx_str[0]  = ' x              '
  amdx_str[1]  = ' y              '
  amdx_str[2]  = ' 1              '
  amdx_str[3]  = ' x^2            '
  amdx_str[4]  = ' x*y            '
  amdx_str[5]  = ' y^2            '
  amdx_str[6]  = ' (x^2+y^2)      '
  amdx_str[7]  = ' x^3            '
  amdx_str[8]  = ' x^2*y          '
  amdx_str[9]  = ' x*y^2          '
  amdx_str[10] = ' y^3            '
  amdx_str[11] = ' x*(x^2+y^2)    '
  amdx_str[12] = ' x*(x^2+y^2)^2  '
                   
  amdy_str[0]  = ' y              '
  amdy_str[1]  = ' x              '
  amdy_str[2]  = ' 1              '
  amdy_str[3]  = ' y^2            '
  amdy_str[4]  = ' y*x            '
  amdy_str[5]  = ' x^2            '
  amdy_str[6]  = ' (x^2+y^2)      '
  amdy_str[7]  = ' y^3            '
  amdy_str[8]  = ' y^2*x          '
  amdy_str[9]  = ' y*x^2          '
  amdy_str[10] = ' x^3            '
  amdy_str[11] = ' y*(x^2+y^2)    '
  amdy_str[12] = ' y*(x^2+y^2)^2  '

  ii = strtrim(indgen(13)+1,2)
  for i=0, 12 do sxaddpar, hdr, 'AMDX'+ii[i], astr.amdx[i], amdx_str[i]
  for i=0, 12 do sxaddpar, hdr, 'AMDY'+ii[i], astr.amdy[i], amdy_str[i]
 
  sxaddhist,'GSSSPUTAST: ' + strmid(systime(),4,20), hdr

  return
end
;------------------------------------------------------------------------------
