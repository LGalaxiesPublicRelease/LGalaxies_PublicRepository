;+
; NAME:
;   dpropmotdisdz
; PURPOSE:
;   Compute differential proper motion distance dD/dz for c=H_0=1.
; CALLING SEQUENCE:
;   dDdz= dpropmotdisdz(z,OmegaM,OmegaL)
; INPUTS:
;   z       - redshift or vector of redshifts
;   OmegaM  - Omega-matter at z=0
;   OmegaL  - Omega-Lambda at z=0
; OPTIONAL INPUTS:
; KEYWORDS
; OUTPUTS:
;   differential proper motion distance in units of the Hubble length c/H_0
; COMMENTS:
; BUGS:
;   May not work for pathological parts of the OmegaM-OmegaL plane.
; EXAMPLES:
; PROCEDURES CALLED:
;   propmotdis
; REVISION HISTORY:
;   2001-Mar-12  Written by Hogg (NYU)
;-
function dpropmotdisdz, z,OmegaM,OmegaL
  TINY= double(1.0e-16)
  ddMdz = 1.0/sqrt((1.0+z)*(1.0+z)*(1.0+OmegaM*z)-z*(2.0+z)*OmegaL)
  OmegaK= 1.0-OmegaM-OmegaL
  if OmegaK LT -TINY then begin
    dM= propmotdis(z,OmegaM,OmegaL)
    ddMdz= sqrt(1.0-OmegaK*dM*dM)*ddMdz
  endif else begin
    if OmegaK GT TINY then begin
      dM= propmotdis(z,OmegaM,OmegaL)
      ddMdz= sqrt(1.0+OmegaK*dM*dM)*ddMdz
    endif
  endelse
  return, ddMdz
end
