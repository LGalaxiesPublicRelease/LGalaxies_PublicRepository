;+
; NAME:
;   near_list
;
; PURPOSE:
;   return index list of positions near a given central position
;
; CALLING SEQUENCE:
;   ind=near_list(racen, deccen, ra, dec, rad, count=count)
;
; INPUTS:
;   racen, deccen - (RA,dec) of central position [deg]
;   ra, dec       - list of coords to test       [deg]
;   radius        - radius to search within      [deg]
;   
; OPTIONAL OUTPUTS:
;   count         - number of elements returned
;
; COMMENTS:
;   Assumes same equinox for input and output coords. 
;   
; REVISION HISTORY:
;   2001-May-25  Written by Douglas Finkbeiner, Princeton
;----------------------------------------------------------------------
function near_list, racen, deccen, ra, dec, rad, count=count

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

; now trim list 
  wdec = where((dec GE dec0) AND (dec LT dec1), dec_ct)
  if dec_ct GT 0 then begin 
     raband = ra[wdec]
     
     IF (ra0 LT ra1) THEN BEGIN 
        wra = where((raband GE ra0) AND (raband LT ra1), ra_ct)
     ENDIF ELSE BEGIN 
        wra = where((raband GE ra0) OR (raband LT ra1), ra_ct)
     ENDELSE 
     
; indices of good vectors
     if ra_ct GT 0 then ind = wdec[wra] else ind = -1
  endif else ind = -1

  count = n_elements(ind)
  if ind[0] EQ -1 then count = 0

  return, ind

end
