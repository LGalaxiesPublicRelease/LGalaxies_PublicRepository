;------------------------------------------------------------------------------
;+
; NAME:
;   djs_posmod
;
; PURPOSE:
;   Return the non-negative modulus x % y, in the range [0,y).
;
; CALLING SEQUENCE:
;   result = djs_posmod(x, y)
;
; INPUTS:
;   xvalue
;
; OUTPUTS:
;   result
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Written D. Schlegel, 15 May 1997, Durham
;-
;------------------------------------------------------------------------------
function djs_posmod, x, y
 
   ; Need 2 parameter
   if N_params() LT 2 then begin
      print, 'Syntax - result = djs_posmod( x, y )'
      return, -1
   endif

;   result = x - y * long(x/y)
   result = x MOD y
   indx = where(result LT 0)
   if (indx[0] NE -1) then result[indx] = result[indx] + abs(y)

   return, result
end
;------------------------------------------------------------------------------
