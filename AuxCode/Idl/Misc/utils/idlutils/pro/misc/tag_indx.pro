;+
; NAME:
;   tag_indx
; PURPOSE:
;   Returns the column number of the given tag to the given structure
;   (-1 if none)
; CALLING SEQUENCE:
;   indx=tag_indx(struct,name)   
; INPUTS:
;   struct - a structure
;   name - name of tag to search for
; OPTIONAL INPUTS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   06-Sep-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function tag_indx, struct, name

return,where(strupcase(tag_names(struct)) eq strupcase(name))

end
