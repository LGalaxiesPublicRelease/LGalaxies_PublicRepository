;+
; NAME:
;   ccm_fitrv.pro
;
; PURPOSE:
;   Does a simple fit to R_V and normalization
;
; CALLING SEQUENCE:
;   R_V = ccm_fitrv(wave, redfac, [weight=, norm=norm, chi2=chi2])
;
; INPUTS:
;   wave   : Wavelength [Ang] at which A is measured
;   redfac : Measured A
;   weight : Weight coefficients (by band)
;
; OUTPUTS:
;   R_V    : Measured R_V
;
; OPTIONAL OUTPUTS:
;   norm   : The measured normalization
;   chi2   : Returned chi2 (unit weights assumed if none passed)
;   
; COMMENTS:
;    The chi2 is measured by func_ccm_fitrv.pro. See that function
;    for weights etc.
;    The reddening curve is computed by ccm_unred.pro
;   
; REVISION HISTORY:
;   2003-Sep-10  Written by D. Finkbeiner & N.Padmanabhan, Princeton
;
;----------------------------------------------------------------------
function ccm_fitrv, wave, redfac, weight=weight, norm=norm, chi2=chi2

  common ccm_fitrv_common, wave_sdss, redfac_meas, red_weight
  wave_sdss   = wave
  redfac_meas = redfac
  if keyword_set(weight) then begin 
     red_weight=weight
  end else red_weight = 0

  par = amoeba(1D-6, function_name='func_ccm_fitrv', P0=[3.1d, 1.0], $
               scale=[2.0d, 1.0d], function_value=chi2vals)

  chi2  = chi2vals[0]
  rv    = par[0]
  norm  = par[1]

  return, rv
end


