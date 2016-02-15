;+
; NAME:
;   is_cap_used
; PURPOSE:
;   Returns whether a cap is used in use_caps
; CALLING SEQUENCE:
;   is_cap_used,use_caps,indx
; INPUTS:
;   use_caps - bit mask indicating which cap is used
;   indx - number indicating which cap we are interestedin
; OPTIONAL INPUTS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
;   Number of caps limited to 32
; PROCEDURES CALLED:
; REVISION HISTORY:
;   09-Nov-2002  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
function is_cap_used,use_caps,indx
return, (use_caps and (2L)^(indx)) gt 0
end
