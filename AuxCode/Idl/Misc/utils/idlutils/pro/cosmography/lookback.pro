;+
; NAME:
;   lookback
; PURPOSE:
;   Compute lookback time (for c/H_0=1).
; CALLING SEQUENCE:
;   t= lookback(z,OmegaM,OmegaL)
; INPUTS:
;   z       - redshift or vector of redshifts
;   OmegaM  - Omega-matter at z=0
;   OmegaL  - Omega-Lambda at z=0
; OPTIONAL INPUTS:
; KEYWORDS
; OUTPUTS:
;   lookback time in units of the Hubble time 1/H_0
; COMMENTS:
; BUGS:
;   The integrator is crude, slow and repetetive.
;   May not work for pathological parts of the OmegaM-OmegaL plane.
; EXAMPLES:
; PROCEDURES CALLED:
;   dlookbackdz()
; REVISION HISTORY:
;   2001-May-12  Written by Hogg (NYU)
;-
;------------------------------------------------------------------------------
function lookback, z,OmegaM,OmegaL

  TINY= double(1.0e-16)
  stepsize= 0.01D              ; minimum stepsize of 0.01
  nsteps= long(z/stepsize)+10  ; minimum of 10 steps
  dz= z/double(nsteps)
  t= double(0.0*z)
  nz= n_elements(z)

  if nz EQ 1 then begin
    if abs(z) LT TINY then t= z else begin
      for zz=0.5*dz,z,dz do t= t+dz*dlookbackdz(zz,OmegaM,OmegaL)
    endelse
  endif else begin
    for i=0L,nz-1L do t[i]= lookback(z[i],OmegaM,OmegaL)
  endelse
  return, t
end
