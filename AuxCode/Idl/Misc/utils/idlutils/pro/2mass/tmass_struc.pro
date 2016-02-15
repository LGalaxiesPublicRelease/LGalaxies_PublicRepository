;+
; NAME:
;   tmass_struc
;
; PURPOSE:
;   Define data structure for FITS copy of 2MASS data
;
; CALLING SEQUENCE:
;   struc = tmass_struc(N)
;
; INPUTS:
;   N     - number of elements in structure array (default 1)
;
; OUTPUTS:
;   struc - 2MASS data structure
;
; COMMENTS:
;   Used by tmass_ascii2fits.  We use DEC instead of the 2MASS
;   name of DECL for declination.
;
; REVISION HISTORY:
;   2003-Jun-26  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function tmass_struc, N

  struc = {tmass_ra:      0.d, $
           tmass_dec:     0.d, $
           tmass_err_maj: 0., $
           tmass_err_min: 0., $
           tmass_err_ang: 0., $
           tmass_j:       0.,$
           tmass_j_ivar:  0., $
           tmass_h:       0.,$
           tmass_h_ivar:  0., $
           tmass_k:       0., $
           tmass_k_ivar:  0., $
           tmass_ph_qual: 'ZZZ', $
           tmass_rd_flg:  0, $
           tmass_bl_flg:  0, $
           tmass_cc_flg:  'zzz', $
           tmass_ndetect: bytarr(3), $
           tmass_nobserve: bytarr(3), $
           tmass_gal_contam: 0b, $
           tmass_mp_flg: 0b, $
           tmass_pts_key: 0L, $
           tmass_hemis: ' ', $
           tmass_jdate: 0d, $
           tmass_dup_src: 0b, $
           tmass_use_src: 0b }

  if keyword_set(N) then struc = replicate(struc, N)

  return, struc
end
