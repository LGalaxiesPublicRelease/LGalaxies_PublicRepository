;+
; NAME:
;   set_use_caps
; PURPOSE:
;   Set the bits in use_caps for a polygon such 
;   that a certain list of caps are being used. Unless
;   /allow_doubles is set, this routine automatically
;   sets use_caps such that no two caps with use_caps
;   set are identical. If /add is set, this routine 
;   doesn't set use_caps to zero before proceeding. 
; CALLING SEQUENCE:
;   set_use_caps,polygon,list [, /allow_doubles, add=add]
; INPUTS:
;   polygon - [Npoly] polygon or polygons to alter
;   list - [Nindices] list of indices to set in each polygon
; OPTIONAL INPUTS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
;  If there are two caps identical except for the sign of cm, we 
;  turn off *all* the caps --- it is zero!
; EXAMPLES:
; BUGS:
;   Number of caps limited to 32
; PROCEDURES CALLED:
; REVISION HISTORY:
;   09-Nov-2002  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
pro set_use_caps, polygon, list, allow_doubles=allow_doubles, add=add, $
                  tol=tol, use_caps=use_caps, $
                  allow_neg_doubles=allow_neg_doubles

if(NOT keyword_set(add)) then polygon.use_caps=0
if(n_elements(tol) eq 0) then tol=1.D-10

for i=0L, n_elements(list)-1L do $
  polygon.use_caps=polygon.use_caps or (2L)^list[i]

if(NOT keyword_set(allow_doubles)) then begin
    for ipoly=0L, n_elements(polygon)-1L do begin
        for i=0L, polygon[ipoly].ncaps-1L do begin
            if(is_cap_used(polygon[ipoly].use_caps,i)) then begin
                for j=i+1L, polygon[ipoly].ncaps-1L do begin
                    if(is_cap_used(polygon[ipoly].use_caps,j)) then begin
                        if(total((*polygon[ipoly].caps)[i].x- $
                                 (*polygon[ipoly].caps)[j].x)^2 lt tol^2) $
                          then begin 
                            if(abs((*polygon[ipoly].caps)[i].cm- $
                                   (*polygon[ipoly].caps)[j].cm) lt tol) then $
                              polygon[ipoly].use_caps= $
                              polygon[ipoly].use_caps and (NOT (2L)^j) $
                            else if(abs((*polygon[ipoly].caps)[i].cm+ $
                                        (*polygon[ipoly].caps)[j].cm) $
                                    lt tol AND $
                                    NOT keyword_set(allow_neg_doubles)) then $
                              polygon[ipoly].use_caps= $
                              polygon[ipoly].use_caps and (NOT (2L)^j) 
                        endif
                    endif
                endfor
            endif
        endfor
    endfor
endif

use_caps=polygon.use_caps

end
