;-----------------------------------------------------------------------
;+
; NAME:
;   djs_angpos
; PURPOSE:
;   Put an angle into the range [0, 360).
;
; CALLING SEQUENCE:
;   result = djs_angpos(xval)
;
; INPUTS:
;   xval
;
; OUTPUTS:
;   result
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Written D. Schlegel, 17 June 1996, Durham
;-
;-----------------------------------------------------------------------
function djs_angpos, xval
 
   retval = (xval MOD 360) + 360 * (xval LT 0)

   return, retval
end
;-----------------------------------------------------------------------
