;+------------------------------------------------------------------------  
; NAME:
;       usno_cone
;+------------------------------------------------------------------------  
; PURPOSE:
;       Determine RA,dec regions to read and call usno_readzone
;+------------------------------------------------------------------------  
; INPUTS:
;   catpath   - catalogue path
;   racen     - RA of region center (J2000)    [degrees]
;   deccen    - dec of region center (J2000)   [degrees]
;   rad       - radius of region               [degrees]
;
; OPTIONAL INPUTS:
;   cattype   - Either 'USNO-A' or 'USNO-B'
;
;+------------------------------------------------------------------------  
; OUTPUTS:
;   result    - float(3,N) array of results.  
;                 data[0,*] = RA (in .01 arcsec)
;                 data[1,*] = (dec+90) (in .01 arcsec)
;                 data[2,*] = magnitudes packed in 32-bit int (see below)
;+------------------------------------------------------------------------  
; COMMENTS:    
;   calls usno_readzone
;+------------------------------------------------------------------------  
; REVISION HISTORY
;   Written  2000 Apr 15 by D. P. Finkbeiner
;   2002-Nov-25  Split from usno_read.pro
;                  and modified to work with USNO-B1.0 - DPF
;+------------------------------------------------------------------------  
;-
PRO usno_cone, catpath, racen, deccen, rad, result, cattype=cattype

  if NOT keyword_set(cattype) then cattype = 'USNO-A'

; -------- For USNO-A (or SA)
  if strupcase(cattype) eq 'USNO-A' then begin 
     rec_len = 12L
     prefix = 'zone'
     swap_if_little_endian = 1B
     zonewidth = 7.5 ; [deg]
  endif 

; -------- For USNO-B
  if strupcase(cattype) eq 'USNO-B' then begin 
     rec_len = 80L
     prefix = 'b'
     swap_if_big_endian = 1B
     zonewidth = 0.1 ; [deg]
     usnob = 1B
  endif 

; declination range
  dec0 = (deccen-rad) > (-90.0)
  dec1 = (deccen+rad) < (+90.0)

; RA range (overly conservative)
  maxdec = abs(dec0) > abs(dec1)
  cosd = cos(maxdec/!radeg)
  IF cosd GT rad/180. THEN BEGIN 
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
  FOR z=z0, z1 DO BEGIN 
     
     zone = long(z*(zonewidth*10))
;     print, z0, z1, zonewidth, zone, z
     if keyword_set(usnob) then begin 
        subdir = string(zone / 10, format='(I3.3)')
     endif else subdir = ''
     path = concat_dir(catpath, subdir)
     IF (ra0 LT ra1) THEN BEGIN 
        usno_readzone, path, zone, ra0, ra1, rec_len, prefix, zdata, $
          swap_if_little_endian=swap_if_little_endian, $
          swap_if_big_endian=swap_if_big_endian
     ENDIF ELSE BEGIN 
        usno_readzone, path, zone, ra0, 360.0, rec_len, prefix, zdata1, $
          swap_if_little_endian=swap_if_little_endian, $
          swap_if_big_endian=swap_if_big_endian
        usno_readzone, path, zone, 0, ra1, rec_len, prefix, zdata2, $
          swap_if_little_endian=swap_if_little_endian, $
          swap_if_big_endian=swap_if_big_endian
        zdata = [[zdata1], [zdata2]]
     ENDELSE 
     
     data = n_elements(data) EQ 0 ? zdata : [[data], [zdata]]
  ENDFOR 

; strip extra dec (works for both USNO-A and B)
  deccat = transpose(data[1, *]) /3.6d5 - 90.
  good = where((deccat GE dec0) AND (deccat LE dec1), ct)
  deccat = 0 ; memory

  IF ct GT 0 THEN BEGIN 
     dtrim = (temporary(data))[*, temporary(good)]
  ENDIF ELSE BEGIN 
     dtrim = temporary(data)  ; keep padding
  ENDELSE 
  
; Now use dot products to strip extras
  racat  = transpose(dtrim[0, *]) /3.6d5
  deccat = transpose(dtrim[1, *]) /3.6d5 - 90.
  uvobj = ll2uv(double([[racat], [deccat]])) ; (n,3) array
  uvcen = ll2uv(double([[racen], [deccen]])) ; (1,3) array
  dot   = temporary(uvobj)#transpose(temporary(uvcen))
  good = where(temporary(dot) GE cos(rad*!dpi/180.d), ct)

  IF ct GT 0 THEN BEGIN 
     result = (temporary(dtrim))[*, good]
  ENDIF ELSE BEGIN 
     result = temporary(dtrim)  ; keep padding
  ENDELSE 

  return
end
