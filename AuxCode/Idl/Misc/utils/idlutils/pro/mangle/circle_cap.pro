;+
; NAME:
;   circle_cap
; PURPOSE:
;   Return a circular cap centered on a certain set of coordinates.
; CALLING SEQUENCE:
;   cap=circle_cap(radius, [ra=, dec=, xyz= ])
; INPUTS:
;   radius - proper radius of cap, in degrees
;   ra, dec - coordinates of center
;    OR
;   xyz - xyz value of center (should be unit vector)
; COMMENTS:
;   Use for the cap of a mangle polygon defining a circle.
;   ra,dec inputs will override xyz inputs
; REVISION HISTORY:
;   01-Oct-2002  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
function circle_cap, ra=ra, dec=dec, xyz=xyz , radius

d2r=!DPI/180.D

if(n_elements(ra) gt 0) then $
  usexyz=angles_to_x(ra,(90.D)-dec) $
else $
  usexyz=xyz

cap=construct_cap()
cap.x=usexyz
cap.cm=(1.D)-cos(radius*d2r)

return,cap

end
