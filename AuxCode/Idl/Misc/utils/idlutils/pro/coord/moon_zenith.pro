;+
; NAME:
;   moon_zenith
;
; PURPOSE:
;   Compute zenith angle of moon, given TAI
;
; CALLING SEQUENCE:
;   zenithang = moon_zenith(TAI, [longitude=, latitude=])
;
; INPUTS:
;   TAI            - time in seconds since MJD 0
;
; OPTIONAL KEYWORDS:
;   longitude      - longitude of observatory [deg] - default to APO
;   latitude       - latitude of observatory [deg]
; 
; OUTPUTS:
;   zenithang      - zenith angle of the moon [deg]
;
;
; COMMENTS:
;   TAI must be specified.
;
; BUGS:
;   Uses geocentric coords, should correct for position on Earth
;
; PROCEDURES CALLED:
;   moonpos
;
; REVISION HISTORY:
;   2001-Apr-06  Written by D. Finkbeiner, APO
;-
;------------------------------------------------------------------------------

function moon_zenith, TAI, longitude=longitude, latitude=latitude

  if n_elements(TAI) EQ 0 then message, 'TAI must be set'
  
; Default to location of Apache Point Observatory
  if (NOT keyword_set(longitude)) then longitude = 360. - 105.820417
  if (NOT keyword_set(latitude)) then latitude = 32.780361

  mjd = TAI/86400.d
  jd = mjd+2400000.5d
  
  MOONPOS, jd, ra, dec

; convert to alt,az
  ct2lst, LST, longitude, junk, jd
  LST = 15. * LST               ; convert from hours to degrees
  HA = LST - ra

  decrad = dec * !dtor
  harad  = HA * !dtor
  latrad = latitude * !dtor

; Compute airmass with spherical trig
  cosz = sin(decrad)*sin(latrad)+ cos(decrad)*cos(harad)*cos(latrad)
  zang = acos(cosz) / !dtor

  return, zang
end
