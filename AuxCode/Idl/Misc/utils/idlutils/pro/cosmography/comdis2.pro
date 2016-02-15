;+
; NAME:
;   comdis2
; PURPOSE:
;   Compute comoving line-of-sight distances (for c/H_0=1).
; CALLING SEQUENCE:
;   D= comdis2(z,OmegaM,OmegaL,cdtable=cdtable,zmaxtable=zmaxtable)
; INPUTS:
;   z       - redshift or vector of redshifts
;   OmegaM  - Omega-matter at z=0
;   OmegaL  - Omega-Lambda at z=0
; OPTIONAL INPUTS:
;   cdtable - returns the lookup table calculated, so you can reuse
;   zmaxtable - set max z value of lookup table
; KEYWORDS
; OUTPUTS:
;   comoving line-of-sight distance in units of the Hubble length c/H_0
; COMMENTS:
;   Call as: 
;      dC = comdis2(z,OmegaM,OmegaL,cdtable=cdtable)
;   if you call it multiple times, so it does not remake the the
;   lookup table.
; BUGS:
;   The integrator is crude, slow and repetitive.
;   May not work for pathological parts of the OmegaM-OmegaL plane.
;   Relies on interpolate() working
; EXAMPLES:
; PROCEDURES CALLED:
;   dcomdisdz()
; REVISION HISTORY:
;   2000-Jun-25  Written by Hogg (IAS)
;   2002-Feb-20  Look-up table to increase speed; Blanton (NYU)
;-
;------------------------------------------------------------------------------
function comdis2, z,OmegaM,OmegaL,cdtable=in_cdtable,zmaxtable=in_zmaxtable

common com_comdis2, cdtable, zmaxtable

if(n_elements(in_cdtable) gt 0) then cdtable=in_cdtable
if(n_elements(in_zmaxtable) gt 0) then zmaxtable=in_zmaxtable

; Set up the table?
  if(NOT keyword_set(zmaxtable)) then zmaxtable=10.d
  resettable=0
  if(n_elements(cdtable) gt 0) then begin
      if(zmaxtable ne cdtable.zmaxtable or $
         OmegaM ne cdtable.OmegaM or $
         OmegaL ne cdtable.OmegaL or $
         max(z) gt cdtable.zmaxtable) then resettable=1
  endif else begin
      resettable=1 
  endelse

; Set it up:
  TINY= double(1.0D-16)
  if(resettable gt 0) then begin
      stepsize= 0.001D           ; minimum stepsize of 0.001
      nsteps= long(zmaxtable/stepsize)+10 ; minimum of 10 steps
      if(zmaxtable lt max(z)) then zmaxtable=1.5D*max(z)
      cdtable={cdstruct, zmaxtable:double(zmaxtable), OmegaM:double(OmegaM), $
               OmegaL:double(OmegaL), nsteps:long(nsteps), $
               zval:dblarr(nsteps), dC:dblarr(nsteps)}
      cdtable.zval=cdtable.zmaxtable*dindgen(nsteps)/double(nsteps)
      hzval=0.5D*(cdtable.zval[1l:nsteps-1l]+cdtable.zval[0l:nsteps-2l])
      dz=cdtable.zval[1]-cdtable.zval[0]
      cdtable.dC[0]=double(0.0D*zmaxtable)
      for iz=1l,cdtable.nsteps-1l do $
          cdtable.dC[iz]= cdtable.dC[iz-1]+ $
        dz*dcomdisdz(hzval[iz-1],OmegaM,OmegaL)
  endif 
      
  dC=interpolate(cdtable.dC,double(cdtable.nsteps)*z/cdtable.zmaxtable)
  return, dC
end
