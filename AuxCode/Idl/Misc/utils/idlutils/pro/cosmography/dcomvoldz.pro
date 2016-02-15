;+
; NAME:
;   dcomvoldz
; PURPOSE:
;   Compute differential comoving volume element dV_c/dz per unit solid angle.
; CALLING SEQUENCE:
;   dVdz= dcomvoldz(z,OmegaM,OmegaL)
; INPUTS:
;   z       - redshift or vector of redshifts
;   OmegaM  - Omega-matter at z=0
;   OmegaL  - Omega-Lambda at z=0
; OPTIONAL INPUTS:
; KEYWORDS
; OUTPUTS:
;   Comoving volume per steradian in units of the Hubble volume (c/H_0)^3
; COMMENTS:
;   Formulae from Carrol, Press & Turner, 1992, Kolb & Turner, 1990, and my
;   own calculation.
; BUGS:
;   May not work for pathological parts of the OmegaM-OmegaL plane.
; EXAMPLES:
; PROCEDURES CALLED:
;   dpropmotdisdz
;   propmotdis
; REVISION HISTORY:
;   2001-Mar-12  Written by Hogg (NYU)
;-
;------------------------------------------------------------------------------
function dcomvoldz, z,OmegaM,OmegaL
  OmegaK= 1.0-OmegaM-OmegaL
  dM= propmotdis(z,OmegaM,OmegaL)
  ddMdz= dpropmotdisdz(z,OmegaM,OmegaL)
  return, (dM*dM*ddMdz/sqrt(1.0+OmegaK*dM*dM))
end
