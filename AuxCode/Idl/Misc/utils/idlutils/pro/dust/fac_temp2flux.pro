;------------------------------------------------------------------------------
;+
; NAME:
;   fac_temp2flux
;
; PURPOSE:
;   Compute factor to convert from brightness temp to flux/sr
;
; CALLING SEQUENCE:
;   fac_temp2flux, nu
;
; INPUTS:
;   nu_ghz -    frequency in GHz
;
; OUTPUTS:
;   <value> -   conversion factor (MJy/sr) / (micro-K)
;
; PROCEDURES CALLED:
;   <none>
;
; COMMENTS:
;
; REVISION HISTORY:
;   Written by D. Finkbeiner, 10 March, 1999 - Berkeley
;   Modified 16 March, 1999 to handle integer nu_ghz argument - DPF
;-
;------------------------------------------------------------------------------

FUNCTION fac_temp2flux, nu_ghz

  k_b=1.3806D-16                ;erg K^-1
  fac = 1d-6*k_b*2/3d10^2*1d18*1d17*double(nu_ghz)^2

  return, fac
END
