;+
; NAME:
;   dlookbackdz
; PURPOSE:
;   Compute differential lookback time dt/dz (for c/H_0=1).
; CALLING SEQUENCE:
;   dtdz= dlookbackdz(z,OmegaM,OmegaL)
; INPUTS:
;   z       - redshift or vector of redshifts
;   OmegaM  - Omega-matter at z=0
;   OmegaL  - Omega-Lambda at z=0
; OPTIONAL INPUTS:
; KEYWORDS
; OUTPUTS:
;   differential lookback time per unit redshift, units of the Hubble time 1/H_0
; COMMENTS:
; BUGS:
;   May not work for pathological parts of the OmegaM-OmegaL plane.
; EXAMPLES:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   2001-May-12  Written by Hogg (NYU)
;-
;------------------------------------------------------------------------------
function dlookbackdz, z,OmegaM,OmegaL
  return, 1.0/((1.0+z)*sqrt((1.0+z)*(1.0+z)*(1.0+OmegaM*z)-z*(2.0+z)*OmegaL)) ;
end
