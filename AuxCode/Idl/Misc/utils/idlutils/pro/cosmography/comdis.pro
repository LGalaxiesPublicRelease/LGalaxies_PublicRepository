;+
; NAME:
;   comdis
; PURPOSE:
;   Compute comoving line-of-sight distances (for c/H_0=1).
; CALLING SEQUENCE:
;   D= comdis(z,OmegaM,OmegaL, weq=weq)
; INPUTS:
;   z       - redshift or vector of redshifts
;   OmegaM  - Omega-matter at z=0
;   OmegaL  - Omega-Lambda at z=0
; OPTIONAL INPUTS:
;   weq     - Equation of state (default=-1)
; KEYWORDS
; OUTPUTS:
;   comoving line-of-sight distance in units of the Hubble length c/H_0
; COMMENTS:
; BUGS:
;   The integrator is crude, slow and repetetive.
;   May not work for pathological parts of the OmegaM-OmegaL plane.
; EXAMPLES:
; PROCEDURES CALLED:
;   dcomdisdz()
; REVISION HISTORY:
;   25-Jun-2000  Written by Hogg (IAS)
;   2004-Sep-06  Added support for different equations of state,
;         Padmanabhan (Princeton) 
;-
;------------------------------------------------------------------------------
function comdis, z,OmegaM,OmegaL, weq=weq

  if (n_elements(weq) EQ 0) then weq=-1.d0

  TINY= double(1.0e-16)
  stepsize= 0.01               ; minimum stepsize of 0.01
  nsteps= long(z/stepsize)+10  ; minimum of 10 steps
  dz= z/double(nsteps)
  dC= double(0.0*z)
  nz= n_elements(z)
  if nz EQ 1 then begin
    if abs(z) LT TINY then dC= z else begin
      for zz=0.5*dz[0],z[0],dz[0] do dC= dC+dz[0]*dcomdisdz(zz,OmegaM,OmegaL, weq=weq)
    endelse
  endif else begin
    for i=0L,nz-1L do dC[i]= comdis(z[i],OmegaM,OmegaL, weq=weq)
  endelse
  return, dC
end
