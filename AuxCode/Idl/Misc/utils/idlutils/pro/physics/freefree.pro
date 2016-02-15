;+
; NAME:
;   freefree
;
; PURPOSE:
;   Calculate freefree spectrum
;
; CALLING SEQUENCE:
;   jfree = freefree(nu, T=T)
;
; INPUTS:
;   nu  - [Hz]
;  
; KEYWORDS:
;   T   - [K] default 10,000K
;   
; OUTPUTS:
;   jfree  - j_nu
;
; OPTIONAL OUTPUTS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   returns answer in erg /cm^3/s/sr/Hz (emissivity, denoted j_\nu)
;   multiply by 1E17 to get MJy/sr/cm 
;      (then times 3E18 to get specific intensity, denoted I_\nu, MJy/sr)
;   multiply by 1E23 to get Jy/sr/cm
;   note for Y=.24 the He correction is about factor 1.42 (for const n_H)
;   
; REVISION HISTORY:
;   2002-Nov-03 - Written by Douglas Finkbeiner, Princeton
;   2006-Dec-01 - moved from Halpha repository - DPF
;
;----------------------------------------------------------------------
function freefree, nu, T=T

  h   = 6.6260755D-27           ; [erg s]
  k_b = 1.380658D-16            ; [erg/K]
  e   = 4.8032068D-10           ; [esu]
  me  = 9.1093897D-28           ; [g]   mass of electron

  Z = 1
  if not keyword_set(T) then T = 1E4
  n_e = 1
  n_i = n_e

  if NOT keyword_set(hefrac) then hefrac = 0.25 ; He mass fraction

; assume HeII / HII ratio is ionrat
  ionrat = 0.02
; -------- He / H number ratio
  herat = hefrac/(4*(1. -hefrac)) * ionrat
  

; -------- H
  n_e = 1. + herat*1  ; electron number density (for n_H = 1)
  n_i = 1             ; H ion number density (for n_H = 1)
  Z = 1.0             ; atomic number of H

  jfree1 = 5.44D-39*gauntff(nu, T)*Z^2*n_e*n_i/T^0.5*exp(-h*nu/k_b/T)

; -------- He
  n_e = 1. + herat*1  ; electron number density (for n_H = 1)
  n_i = herat         ; He ion number density (for n_H = 1)
  Z = 1.0             ; charge of Helium (singly ionized)

  jfree2 = 5.44D-39*gauntff(nu, T)*Z^2*n_e*n_i/T^0.5*exp(-h*nu/k_b/T)
  jfree = jfree1+jfree2

  return, jfree
end


pro rayleigh2em
  alpha32 = 11.7D-14            ; Spitzer; 10000K

  h   = 6.6260755D-27           ; erg s
  nu = 3D18/6563.               ; Hz

  j = h*nu*alpha32/4/!dpi       ; erg/cm^3/s/sr

  pc = 86400.d*365.2425 * 29979245800.D * 3.2616
  R = alpha32 * pc * 1E-6       ; 1E6/4!pi phot/cm^2/s/sr (1EM in R)

  return
end

