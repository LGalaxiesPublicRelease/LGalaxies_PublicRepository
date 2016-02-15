;+
; NAME:
;       tmass_read
;
; PURPOSE:
;       Read 2MASS stars near a given RA, dec
;
; INPUTS:
;   racen     - RA of region center (J2000)    [degrees]
;   deccen    - dec of region center (J2000)   [degrees]
;   rad       - radius of region               [degrees]
;
;
; OUTPUTS:
;   result    - data structure defined by FITS file
;
; COMMENTS:    
;   Determines RA,dec regions to read and calls tmass_readzone
;
; REVISION HISTORY
;   2003-Jul-14   Written by D. P. Finkbeiner
;   2005-Aug-18   Return 0 if no objects are within radius - DPF
;
;------------------------------------------------------------------------  
function tmass_read, racen, deccen, rad

  fitspath = concat_dir(getenv('TWOMASS_DIR'), 'slice')
  if fitspath EQ 'slice' then begin 
     splog, 'You must set your $TWOMASS_DIR environment variable!'
     stop
  endif

  prefix = '2mass-'
  zonewidth = 0.1               ; [deg]

; declination range
  dec0 = (deccen-rad) > (-90.0)
  dec1 = (deccen+rad) < (+90.0)

; RA range (overly conservative)
  maxdec = abs(dec0) > abs(dec1)
  cosd = cos(maxdec/!radeg)
  IF cosd GT rad/180. THEN BEGIN ; avoid erroneous wrap
     ra0 = ((racen-rad/cosd) + 360) MOD 360
     ra1 = ((racen+rad/cosd) + 360) MOD 360
  ENDIF ELSE BEGIN 
     ra0 = 0.
     ra1 = 360.
  ENDELSE
     
; dec zone numbers
  z0 = floor((90+dec0)/zonewidth) < long(180.0/zonewidth - 1)
  z1 = floor((90+dec1)/zonewidth) < long(180.0/zonewidth - 1)

; loop over zones
  FOR zone=z0, z1 DO BEGIN 
     
     subdir = string(zone / 10, format='(I3.3)')
     path = concat_dir(fitspath, subdir)
     IF (ra0 LT ra1) THEN BEGIN 
        tmass_readzone, path, zone, ra0, ra1, prefix, zdata
     ENDIF ELSE BEGIN 
        tmass_readzone, path, zone, ra0, 360.0, prefix, zdata1
        tmass_readzone, path, zone, 0, ra1, prefix, zdata2
        zdata = [zdata1, zdata2]
     ENDELSE 
     
     data = n_elements(data) EQ 0 ? zdata : [data, zdata]
  ENDFOR 

; strip extra dec (works for both TMASS-A and B)
  racat  = data.tmass_ra
  deccat = data.tmass_dec
  good = where((deccat GE dec0) AND (deccat LE dec1), ct)

  IF ct GT 0 THEN BEGIN 
     dtrim = data[good]
  ENDIF ELSE BEGIN 
     return, 0B         ; give up and return nothing
  ENDELSE 
  
; Now use dot products to strip extras
  racat  = dtrim.tmass_ra
  deccat = dtrim.tmass_dec
  uvobj  = ll2uv(double([[racat], [deccat]]), /double) ; (n,3) array
  uvcen  = ll2uv(double([[racen], [deccen]]), /double) ; (1,3) array
  dot    = uvobj#transpose(uvcen)
  good   = where(dot GE cos(rad*!dpi/180.d), ct)

  IF ct GT 0 THEN BEGIN 
     result = dtrim[good]
  ENDIF ELSE BEGIN 
     return, 0B         ; give up and return nothing
  ENDELSE 

  return, result
end
