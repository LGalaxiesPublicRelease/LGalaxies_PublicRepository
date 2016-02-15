;+
; NAME:
;   lumdis
; PURPOSE:
;   Compute luminosity distances (for c/H_0=1).
; CALLING SEQUENCE:
;   D= lumdis(z,OmegaM,OmegaL, weq=weq)
; INPUTS:
;   z       - redshift or vector of redshifts
;   OmegaM  - Omega-matter at z=0
;   OmegaL  - Omega-Lambda at z=0
; OPTIONAL INPUTS:
;   weq     - Equation of state (default=-1)
; KEYWORDS
; OUTPUTS:
;   luminosity distance in units of the Hubble length c/H_0
; COMMENTS:
; BUGS:
;   May not work for pathological parts of the OmegaM-OmegaL plane.
; EXAMPLES:
; PROCEDURES CALLED:
;   propmotdis()
; REVISION HISTORY:
;   25-Jun-2000  Written by Hogg (IAS)
;   2004-Sep-8, Added equation of state for OmegaL, Padmanabhan
;   (Princeton)
;-
function lumdis, z,OmegaM,OmegaL, weq=weq
  return, propmotdis(z,OmegaM,OmegaL,weq=weq)*(1.0+z)
end
