;------------------------------------------------------------------------------
;+
; NAME:
;   fac_flux2temp
;
; PURPOSE:
;   Compute factor to convert from flux/sr to brightness temp
;
; CALLING SEQUENCE:
;   fac_flux2temp, nu_ghz
;
; INPUTS:
;   nu_ghz -    frequency in GHz
;
; OUTPUTS:
;   <value> -   conversion factor (micro-K) / (MJy/sr)
;
; PROCEDURES CALLED:
;   fac_temp2flux()
;
; COMMENTS:
;   see fac_temp2flux.pro 
;   We call fac_temp2flux so that these routines are the inverse of
;   each other by construction. 
;
; REVISION HISTORY:
;   Written by D. Finkbeiner, 10 March, 1999 - Berkeley
;-
;------------------------------------------------------------------------------

FUNCTION fac_flux2temp, nu_ghz

  return, 1./fac_temp2flux(nu_ghz)
END 
