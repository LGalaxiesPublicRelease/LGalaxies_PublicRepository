;+
; NAME:
;   legendre_zero
;
; PURPOSE:
;   Compute zeros of the m=0 Legendre polynomials
;
; CALLING SEQUENCE:
;   thzero = legendre_zero(l)
;
; INPUTS:
;   l - order of Legendre polynomial Pl0(x) to compute
;   
; KEYWORDS:
;   
; OUTPUTS:
;   thzero  - theta values of zeros [radians]
;
; EXAMPLES:
;   
; COMMENTS:
;   Just calls fx_root.  We need the common block and the function
;     legfn. 
;
; REVISION HISTORY:
;   2003-Feb-21  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function legendre_zero, l

; -------- store the degree of the legendre polynomial in common block
  common legendre_degree, deg

  deg = l

; -------- equally spaced guesses in theta
  thguess = (dindgen(l)+0.5)/l*!dpi
  nzero = l

  worst = 0
  thzero = dblarr(nzero)
  for i=0L, nzero/2 do begin 
     guess = [-1, 0, 1]*!dpi/l/4+thguess[i]
     x = fx_root(cos(guess), 'legendre_zero_func', /double, tol=1d-12, itmax=1000)

     check = legendre(x, l, 0, /double)
     worst = worst > abs(check)
     thzero[i] = acos(x)
  endfor 
  if nzero GT 1 then $
    thzero[nzero-nzero/2:nzero-1] = !dpi-reverse(thzero[0:nzero/2-1])
;  print, 'Worst case', worst

  return, thzero
end
