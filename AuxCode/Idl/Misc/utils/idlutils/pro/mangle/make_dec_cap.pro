;+
; NAME:
;   make_dec_cap
; PURPOSE:
;   Creates a structure containing a cap expressing a limit in dec
; CALLING SEQUENCE:
;   cap=make_dec_cap(declimit, [,sign= ]
; INPUTS:
;   declimit - limit on dec 
; OPTIONAL INPUTS:
;   sign - sign of the cap (default 1.)
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
function make_dec_cap,declimit,sign=sign

if(n_params() NE 1) then $
  message, 'Usage: cap=make_dec_cap(declimit [,sign= ]) '

if(n_elements(sign) eq 0) then sign=1.D

; make cap structure
cap1=construct_cap()
cap=replicate(cap1,n_elements(declimit))

; get direction of cap
ra=0.
dec=90.
cap.x=angles_to_x(ra,(90.D)-dec)

; set size of cap
cap.cm=(1.-sin(declimit*!dpi/180.D))*((sign gt 0.) ? 1. : -1.)

return,cap

end
