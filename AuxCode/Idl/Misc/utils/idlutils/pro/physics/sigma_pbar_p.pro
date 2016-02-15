;+
; NAME:
;   sigma_pbar_p
;
; PURPOSE:
;   Annihilation cross section of pbar on p as a function of sqrt(s)
;   
; CALLING SEQUENCE:
;   sigma = sigma_pbar_p(E_in)
;
; INPUTS:
;   E_in    - Input energy, default is sqrt(s) [MeV]
;             if /labframe set, then total energy of pbar in p rest-frame.
;
; KEYWORDS:
;   labframe - set to interpret energy argument as total energy of
;              pbar in rest from of proton, appropriate for
;              anti-proton cosmic rays.
;
; OUTPUTS:
;   sigma    - cross section in barns.  (1 barn = 1E-24 cm^2)
;
; OPTIONAL OUTPUTS:
;   
; RESTRICTIONS:
;   If total energy is LESS than the rest mass of the two particles,
;   print warnings.
;
; EXAMPLES:
;     sig = sigma_pbar_p(sqrts)
;
; COMMENTS:
;   Cross section formula taken from
;     http://atlas.web.cern.ch/Atlas/GROUPS/SOFTWARE/OO/domains/simulation/G4PhysicsStudies/documentation/ameline/node70.html
;   sqrt(s) is in MeV and corresponds to center of mass total energy
;   (including rest mass of the two particles)
;   s = (p1+p2)^2 where p1 and p2 are the 4-vectors (E,px,py,pz), and 
;    the "square" operation uses the Minkowski metric (1,-1,-1,-1) and c=1.
;   
; REVISION HISTORY:
;   2006-Feb-04  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function sigma_pbar_p, E_in, labframe=labframe

  mN = 939.566   ; [MeV] mass of neutron
  sigma0 = 0.120 ; [barn]
  s0 = 4*mN^2    ; 4 times mass of neutron squared
  A = 50.        ; [MeV]
  B = 0.6        ; dimensionless parameter

  if keyword_set(labframe) then begin
     s = 2*mN*(mN+E_in)     ; assume E_in is labframe energy
  endif else begin 
     s = E_in*E_in          ; assume E_in is sqrt(s)=c.o.m energy
  endelse
  w0 = where((s LT s0) OR (E_in LE 0), n0)
  if n0 GT 0 then begin
     splog, 'energy underflow - filling in with zero cross-section'
     s[w0] = s0  ; fill in with harmless value
  endif

  sigma = sigma0*(s0/s)*(A*A*s0/((s-s0)^2 + A*A*s0) + B)

  if n0 GT 0 then begin
     sigma[w0] = 0   ; zero out array elements for bad input values
  endif

  return, sigma
end


pro test

  m = 939.566 ; MeV
  Ej = findgen(10000)+m
  sqrts = sqrt(2*m*(m+Ej))
  sig = sigma_pbar_p(sqrts)

  plot, Ej, sig

  return
end
