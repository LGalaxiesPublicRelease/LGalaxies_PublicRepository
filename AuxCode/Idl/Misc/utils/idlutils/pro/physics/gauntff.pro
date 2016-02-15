;+
; NAME:
;   gauntff
;
; PURPOSE:
;   Compute the Gaunt factor for free-free
;
; CALLING SEQUENCE:
;   gff = gauntff(nu, T)
;
; INPUTS:
;   nu       - frequency   [Hz]
;   T        - temperature [Kelvin]
;
; EXAMPLES:
;   gff = gauntff(1.0D10, 1D4)
;   
; COMMENTS:
;   This approximation is valid for frequencies far above the plasma
;    frequency, but far below kT/h. 
;
; REVISION HISTORY:
;   2002-Jun-14   Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function gauntff, nu, T

; -------- The constants, from http://physics.nist.gov, 1998 CODATA
  h   = 6.62606876D-27          ; erg s   (+/- 0.000 000 52 D-27)
  k_b = 1.3806503D-16           ; erg/K   (+/- 0.000 0024 D-16)
  c   = 299792458.D             ; m/s     (exact)
  e_C = 1.602176462D-19         ; Coulomb (+/- 0.000 000 063 D-19)
  e   = e_C * c *10             ; esu     (see Jackson)

  mp   = 1.67262158D-24         ; g  (proton mass) (+/- 0.000 000 13 D-24)
  me   = 9.10938188D-28         ; g  (electron mass) (+/- 0.000 000 72 D-28)
  euler = 0.57721566490153286D  ; don't say that isn't enough precision!
  Z = 1


; -------- The Gaunt factor from Spitzer, p. 58 (eq 3-55). 
  gff = sqrt(3.d)/!dpi * ( $
    0.5*alog(8*k_b^3*double(T)^3/(!dpi^2 * Z^2 * e^4 * me * double(nu)^2)) - $
                           2.5d * euler)

  return, gff
end
