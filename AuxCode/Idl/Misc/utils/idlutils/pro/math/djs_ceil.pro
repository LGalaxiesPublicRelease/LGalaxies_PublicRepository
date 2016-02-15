;------------------------------------------------------------------------------
;+
; NAME:
;   djs_ceil
;
; PURPOSE:
;   Return smallest integer not less than xvalue.
;   This is identical to the C library function "ceil()".
;
; CALLING SEQUENCE:
;   result = djs_ceil(xvalue)
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
;   Written D. Schlegel, 27 November 1996, Durham
;-
;------------------------------------------------------------------------------
function djs_ceil, xvalue
 
   ; Need 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - result = djs_ceil( xvalue )'
      return, -1
   endif

   result = - long(-xvalue + long(xvalue) + 1) + long(xvalue) + 1

   return, result
end
;------------------------------------------------------------------------------
