;+
; NAME:
;   make_ra_cap
; PURPOSE:
;   Creates a structure containing a cap expressing a limit in ra
; CALLING SEQUENCE:
;   cap=make_ra_cap(ralimit, [,sign= ])
; INPUTS:
;   ralimit - limit on ra 
; OPTIONAL INPUTS:
;   sign - sign of the cap (default 1.)
; REVISION HISTORY:
;   01-Oct-2002  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
function make_ra_cap,ralimit,sign=sign

if(n_params() NE 1 ) then $
  message, 'Usage: cap=make_ra_cap(ralimit [,sign=, node=,]) '

if(n_elements(sign) eq 0) then sign=1.D

; make cap structure
cap1=construct_cap()
cap=replicate(cap1,n_elements(ralimit))

; get direction of cap
ra=ralimit+90. mod 360.
dec=0.+dblarr(n_elements(ralimit))
cap.x=angles_to_x(ra,(90.D)-dec)

; set size of cap
cap.cm=1.*((sign gt 0.) ? 1. : -1.)

return,cap

end
