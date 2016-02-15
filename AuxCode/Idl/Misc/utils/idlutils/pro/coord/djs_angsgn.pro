;-----------------------------------------------------------------------
;+
; NAME:
;   djs_angsgn
; PURPOSE:
;  Put an angle into the range [-180, 180).
;
; CALLING SEQUENCE:
;   result = djs_angsgn(xval)
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
function djs_angsgn, xval

   retval = (xval MOD 360.d0) + 360.d0 * (xval LT 0.d0)
   retval = retval - 360.d0 * (retval GE 180.d0)

   return, retval
end
;-----------------------------------------------------------------------
