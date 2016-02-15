;+
; NAME:
;   copy_caps
; PURPOSE:
;   copy information about caps from one polygon to another
; CALLING SEQUENCE:
;   copy_caps, poly1, poly2
; INPUTS:
; OPTIONAL INPUTS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
;   blows away old pointers
; EXAMPLES:
; BUGS:
;   Number of caps limited to 32
; PROCEDURES CALLED:
; REVISION HISTORY:
;   01-Oct-2002  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
pro copy_caps,poly1,poly2

for i=0L, n_elements(poly1)-1L do begin
    poly2[i].ncaps=poly1[i].ncaps
    poly2[i].use_caps=poly1[i].use_caps
    sz=n_elements(*poly1[i].caps)
    tmp=replicate(construct_cap(),sz)
    struct_assign,(*poly1[i].caps),tmp
    poly2[i].caps=ptr_new(tmp)
endfor

return
end
