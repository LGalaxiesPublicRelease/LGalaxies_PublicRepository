;+
; NAME:
;   dcomdisdz
; PURPOSE:
;   Compute differential comoving line-of-sight distances (for c/H_0=1).
; CALLING SEQUENCE:
;   dDdz= dcomdisdz(z,OmegaM,OmegaL, weq=weq)
; INPUTS:
;   z       - redshift or vector of redshifts
;   OmegaM  - Omega-matter at z=0
;   OmegaL  - Omega-Lambda at z=0
; OPTIONAL INPUTS:
;   weq     - Eq. of state (Default =-1)
; KEYWORDS
; OUTPUTS:
;   differential comoving distance DD/dz in units of the Hubble length c/H_0
; COMMENTS:
; BUGS:
;   May not work for pathological parts of the OmegaM-OmegaL plane.
; EXAMPLES:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   25-Jun-2000  Written by Hogg (IAS)
;   2004-Sep-06  Added support for different equations of state,
;         Padmanabhan (Princeton)
;-
function dcomdisdz, z,OmegaM,OmegaL, weq=weq

  if (n_elements(weq) EQ 0) then weq=-1.d0

  inva = (1.d0 + z)
  return, 1.d0/sqrt(inva*inva*(OmegaM*z+1.d0-OmegaL) + OmegaL*(inva^(3.0d0*(1.d0+weq))))

end
