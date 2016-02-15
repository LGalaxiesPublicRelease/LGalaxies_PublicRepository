;+
; NAME:
;   construct_cap
; PURPOSE:
;   Create the structure for a cap
; CALLING SEQUENCE:
;   poly=construct_cap()
; INPUTS:
; OPTIONAL INPUTS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   01-Oct-2002  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
function construct_cap

; make cap structure
cap1={capstr, x:dblarr(3), cm:0.D}

return,cap1

end
