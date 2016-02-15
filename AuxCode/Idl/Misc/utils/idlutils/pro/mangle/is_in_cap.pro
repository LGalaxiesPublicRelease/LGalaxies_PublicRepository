;+
; NAME:
;   is_in_cap
; PURPOSE:
;   Is an xyz (or radec) position in a given cap?
; CALLING SEQUENCE:
;   result=is_in_cap(ra=, dec=, xyz=, cap )
; INPUTS:
;   ra - set of ra values
;   dec - set of dec values
;   xyz - xyz value(s) (overrides ra and dec)
;   cap - single cap to check
; OPTIONAL INPUTS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
;   Either ra and dec, or xyz must be set; xyz overrides ra and dec
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   01-Oct-2002  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
function is_in_cap, ra=ra, dec=dec, xyz=xyz, cap

if(n_elements(xyz) gt 0) then $
  usexyz=xyz $
else $
  usexyz=angles_to_x(ra,(90.D)-dec)

nxyz=n_elements(usexyz)/3L
dotp=(transpose(reform(usexyz,3,nxyz))#reform(cap.x,3,1))
if(cap.cm gt 0) then $
  return,1.-dotp lt abs(cap.cm) $
else $
  return,1.-dotp gt abs(cap.cm) 

end
