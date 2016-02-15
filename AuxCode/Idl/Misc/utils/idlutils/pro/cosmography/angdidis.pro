;+
; NAME:
;   angdidis
; PURPOSE:
;   Compute angular diameter distancea (for c/H_0=1).
; CALLING SEQUENCE:
;   D= angdidis(z,OmegaM,OmegaL, weq=weq)
; INPUTS:
;   z       - redshift or vector of redshifts
;   OmegaM  - Omega-matter at z=0
;   OmegaL  - Omega-Lambda at z=0
; OPTIONAL INPUTS:
;   weq     - Equation of state (Default=-1)
; KEYWORDS
; OUTPUTS:
;   angular diameter distance in units of the Hubble length c/H_0
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
function angdidis, z,OmegaM,OmegaL, weq=weq
  return, propmotdis(z,OmegaM,OmegaL,weq=weq)/(1.0+z)
end
