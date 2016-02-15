;+
; NAME:
;   evolve_electron_spectrum
;
; PURPOSE:
;   Evolve an initial electron spectrum to a later time
;
; CALLING SEQUENCE:
;   spec = evolve_electron_spectrum(logE, initspec, logEout, time)
;
; INPUTS:
;   logE     - log10 input energy bins [GeV]
;   initspec - initial spectrum, E^2 dN/dE [erg]
;   logEout  - log10 output energy bins [GeV]
;   time     - time since t0 [years]
;
; OUTPUTS:
;   spec     - evolved spectrum, E^2 dN/dE [erg]
;
; EXAMPLES:
;   see evolve_powerlaw_electrons.pro
;
; COMMENTS:
;   Energy is log10 but spectrum is linear!!!
;
; REVISION HISTORY:
;   2007-Feb-06 - Taken from grb_electron_spec.pro 
;                           by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
function evolve_electron_spectrum, logE, initspec, logEout, time

; -------- set defaults
  if not keyword_set(logE)     then stop
  if not keyword_set(initspec) then stop
  if not keyword_set(logEout)  then logEout = logE
  if not keyword_set(time)     then time = 0.d

; -------- set parameters
  gamma = (10.d ^logEout) / 511.d-6 ; gamma of each final energy bin
  dlngamma  = alog(10) * (logEout[1]-logEout[0]) ; d ln(gamma)

  tau = time/1d11               ; [time is in yr]

  gamma0 = gamma/(1 - gamma*tau) ; electrons with gamma now had gamma0
  w = where((gamma*tau) GE 1, nw)
  if nw GE 1 then gamma0[w] = 1d100  ; ???
  print, 'nw', nw

  dlngamma0 = dlngamma*(gamma0/gamma)

; -------- interpolate off of initspec
  logE0 = alog10(gamma0/gamma)+logEout
  thisspec = interpol(alog10(double(initspec) > 1D-300), logE, logE0)

  erg0 = 10.d ^ logE0 / 624.150974d ; original binning energy [erg]
  erg = 10.d ^ logEout / 624.150974d   ; result binning energy [erg]
  NdlnE = (10.d ^ thisspec) /erg0 
  
  N1dlnE = NdlnE*dlngamma0/dlngamma
  N1 = erg*N1dlnE                   ; E^2 dN/dE for the result
  if nw GE 1 then n1[w] = 1d-300

  return, N1
end
