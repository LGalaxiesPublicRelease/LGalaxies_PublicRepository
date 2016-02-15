;+
; NAME:
;   func_ccm_fitrv, params
;
; PURPOSE:
;   compute chi2 for the ccm_fitrv function
;
; CALLING SEQUENCE:
;   called via a fitting program -- see ccm_fitrv.pro
;
; INPUTS:
;   params   -   R_V  = params[0];   norm = params[1]
;
; OUTPUTS:
;   chi2     -  chi squared for fit. 
;
; COMMENTS:
;   wavelength, redfac, and weights are passed via common block. 
;   unit weights are assumed if none are passed. 
;
; REVISION HISTORY:
;   2003-Sep-10  Written by D. Finkbeiner & N.Padmanabhan, Princeton
;
;----------------------------------------------------------------------
function func_ccm_fitrv, params

  if n_elements(params) ne 2 then stop
  R_V  = params[0]
  norm = params[1]

  common ccm_fitrv_common, wave_sdss, redfac_meas, red_weight
  
  if NOT keyword_set(red_weight) then $
    red_weight = 1+fltarr(n_elements(wave_sdss)) 

  if n_elements(wave_sdss) NE n_elements(redfac_meas) then stop

  ccm_unred, wave_sdss, fltarr(5)+1, 1, funred, R_V=R_V
; The 1.0519 comes from the SFD definition of B and V (and Ebv)
  red_R_V = 2.5*alog10(funred)*1.0519
  
  chi2 = total((red_R_V-redfac_meas*norm)^2*red_weight)/total(red_weight)

  return, chi2
end
