function asinh2, x
;+
; NAME:
;     ASINH
; PURPOSE:
;     Return the inverse hyperbolic sine of the argument
; EXPLANATION:
;     The inverse hyperbolic sine is used for the calculation of asinh 
;     magnitudes, see Lupton et al. (1999, AJ, 118, 1406)
;
; CALLING SEQUENCE
;     result = asinh( x) 
; INPUTS:
;     X - hyperbolic sine, numeric scalar or vector (not complex) 
;
; OUTPUT:
;     result - inverse hyperbolic sine, same number of elements as X
;              double precision if X is double, otherwise floating pt.
;
; METHOD:
;     Expression given in  Numerical Recipes, Press et al. (1992), eq. 5.6.7 
;     Note that asinh(-x) = -asinh(x) and that asinh(0) = 0. and that
;     if y = asinh(x) then x = sinh(y).     
;
; REVISION HISTORY:
;     Written W. Landsman                 February, 2001
;     Bug fixed for multiple elements  M Blanton Mar 2002
;-

 y = x*0.      ;Output vector is at least floating pt.
 y = alog( abs(x) + sqrt( x^2 + 1.d0) )

 index = where(x LT 0 ,count)
 if count GT 0 then y[index] = -y[index]

 return, y

 end
