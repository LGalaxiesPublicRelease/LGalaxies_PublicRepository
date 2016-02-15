;+
; NAME:
;   dpropdisdz
; PURPOSE:
;   Compute differential proper line-of-sight distances (for c/H_0=1).
; CALLING SEQUENCE:
;   dDdz= dpropdisdz(z,OmegaM,OmegaL)
; INPUTS:
;   z       - redshift or vector of redshifts
;   OmegaM  - Omega-matter at z=0
;   OmegaL  - Omega-Lambda at z=0
; OPTIONAL INPUTS:
; KEYWORDS
; OUTPUTS:
;   differential proper distance DD/dz in units of the Hubble length c/H_0
; COMMENTS:
; BUGS:
;   May not work for pathological parts of the OmegaM-OmegaL plane.
; EXAMPLES:
; PROCEDURES CALLED:
;   dcomdisdz()
; REVISION HISTORY:
;   25-Jun-2000  Written by Hogg (IAS)
;-
function dpropdisdz, z,OmegaM,OmegaL
  return, (dcomdisdz(z,OmegaM,OmegaL)/(1.0+z))
end
