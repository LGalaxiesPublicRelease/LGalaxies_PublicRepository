;+
; NAME:
;   lg_to_helio
; PURPOSE:
;   Convert to heliocentric redshifts from Local-Group-centric redshifts.
; CALLING SEQUENCE:
;   z_helio = lg_to_helio(z_lg,ra,dec)
; INPUTS:
;   z_lg   - local-group-centric redshift
;   RA        - right ascension (deg, J2000)
;   Dec       - declination (deg, J2000)
; OUTPUTS:
;   z_helio      - heliocentric redshift
; REVISION HISTORY:
;   Originally imported by Hogg, 2002-08 or so.
;   MRB corrected sign error in correction 2004-04-08 (affected tags
;   v5_0_0 and previous)
;-
function lg_to_helio, z_lg,RA,Dec

; from Yahil, 1977, ApJ, 217, 903
  z_sun_lg = 308.0D / 300000.0D
  RA_lg = 343D
  Dec_lg = 52D

  z_h = z_lg - z_sun_lg * cos(double(!PI) * djs_diff_angle(RA,Dec,RA_lg,Dec_lg) / 180D)
  return, z_h
end
