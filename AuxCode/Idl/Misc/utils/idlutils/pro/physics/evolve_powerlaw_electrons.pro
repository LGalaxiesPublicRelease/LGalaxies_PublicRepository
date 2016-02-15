;+
; NAME:
;   evolve_powerlaw_electrons
;
; PURPOSE:
;   Evolve a power-law CR electron spectrum to a later time
;
; CALLING SEQUENCE:
;   spec = evolve_powerlaw_electrons(ind, time=time, Ebin=Ebin)
;
; INPUTS:
;   ind   - power law index
;   time  - time [yr] after injection at which to evaluate spectrum
;
; OUTPUTS:
;   Ebin  - Energy bins [GeV]
;   spec  - the evolved electron spectrum, in E^2 dN/dE [erg]
;   
; EXAMPLES:
;   
; COMMENTS:
;   
; REVISION HISTORY:
;   2007-Feb-06  Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
function evolve_powerlaw_electrons, ind, time=time, Ebin=Ebin

  logE     = (dindgen(600)/100)-1.25 ; [log GeV]
  Ebin     = 10.d ^ logE             ; [GeV]
  initspec = Ebin^ind                ; initial spectrum, E^2 dN/dE [erg]

  spec = evolve_electron_spectrum(logE, initspec, logE, time)

  return, spec
end
