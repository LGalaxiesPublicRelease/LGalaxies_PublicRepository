;+
; NAME:
;   cap_distance
;
; PURPOSE:
;   Return distance from coordinates to a cap, in degrees.
;
; CALLING SEQUENCE:
;   cap_distance
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The sign is positive if in the cap, and negative if outside
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   angles_to_x()
;
; REVISION HISTORY:
;   19-Jun-2003  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function cap_distance, cap, xyz=xyz, ra=ra, dec=dec

if(n_tags(cap) eq 0 or $
   (n_elements(xyz) eq 0 and $
    (n_elements(ra) eq 0 or n_elements(dec) eq 0))) then begin
    print, 'Usage: dist= cap_distance(cap, [xyz=, ra=, dec=]) '
    print, '  (where cap is a CAP structure and either xyz OR '
    print, '   ra and dec are set)'
    return,-9999.
endif

if(n_elements(xyz) eq 0) then $
  xyz = angles_to_x(ra,(90.D)-dec)

nobj = n_elements(xyz) / 3L
dotprod = transpose(reform(xyz,3,nobj)) # reform(cap.x,3,1)
cdist = (acos((1.D) - abs(cap.cm)) - acos(dotprod)) * 180.D0 / !dpi
cdist = cdist - 2 * cdist * (cap.cm LT 0) ; Flip the sign if CM is negative

   return, cdist
end
