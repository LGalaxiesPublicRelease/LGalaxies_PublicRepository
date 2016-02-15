; NAME:
;   legendre_zero_func
;
; PURPOSE:
;   Function for fx_root to call from legendre_zero.pro
;
; REVISION HISTORY:
;   2003-Feb-21  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------

function legendre_zero_func, x
  common legendre_degree, deg
  return, legendre(x, deg, 0, /double)
end
