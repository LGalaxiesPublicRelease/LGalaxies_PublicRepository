;------------------------------------------------------------------------------
;+
; NAME:
;   djs_mean
;
; PURPOSE:
;   Return the mean value of an array.
;
; CALLING SEQUENCE:
;   result = djs_mean(array)
;
; INPUTS:
;   array      - Array of numbers
;
; OUTPUTS:
;   result     - Computed mean
;
; PROCEDURES CALLED:
;
; COMMENTS:
;   This function is faster than the IDL function MEAN(), and will not crash
;   when passed a 1-element array.
;
; REVISION HISTORY:
;   06-Oct-1997  Written by David Schlegel, Durham.
;-
;------------------------------------------------------------------------------
function djs_mean, array
 
   ; Need 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - result = djs_mean( array )'
      return, -1
   endif

   result = total( array + 0D ) / N_elements(array)

   return, result
end
;------------------------------------------------------------------------------
