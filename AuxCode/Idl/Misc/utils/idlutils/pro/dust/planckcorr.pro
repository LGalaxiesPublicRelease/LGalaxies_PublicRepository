;------------------------------------------------------------------------------
;+
; NAME:
;   planckcorr
; PURPOSE:
;   Compute factor to convert from brightness temp to thermodynamic temp
;
; CALLING SEQUENCE:
;   planckcorr, nu_ghz
;
; INPUTS:
;   nu_ghz -    frequency in GHz
;
; OUTPUTS:
;   <value> -   conversion factor
;
; PROCEDURES CALLED:
;   <none>
;
; COMMENTS:
;   This conversion factor corresponds to the PLNCKCOR FITS header
;   keyword found in COBE/DMR data files produced by NASA.
;   For comparison, their results for 31.5, 53, and 90 GHz are
;
;   PLNCKCOR=             1.025724 /  Thermodynamic temperature = 
;   PLNCKCOR=             1.074197 /  Thermodynamic temperature =  
;   PLNCKCOR=             1.225941 /  Thermodynamic temperature =
;   COMMENT                        /   PLNCKCOR * antenna temperature  
;
; REVISION HISTORY:
;   Written by D. Finkbeiner, 10 March, 1999 - Berkeley
;-
;------------------------------------------------------------------------------

FUNCTION planckcorr, nu_ghz

  k_b = 1.3806E-23              ; J/K
  h = 6.6262E-34                ; J*s
  T_cmb = 2.73                  ; K	
  
  nu = nu_ghz*1.0E9             ; Hz
  x = h * nu / (k_b * T_cmb)
  
  result = (exp(x)-1.)^2. / (x^2. * exp(x)) ; thermodynamic uK  

  return, result
END
