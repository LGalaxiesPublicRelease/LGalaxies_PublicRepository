;+
; NAME:
;   destruct_polygon
; PURPOSE:
;   Destroy the structure for a polygon. 
; CALLING SEQUENCE:
;   destruct_polygon, poly
; INPUTS:
;   poly - polygon to destroy
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
pro destruct_polygon,polygon

for i=0L, n_elements(polygon)-1L do $
  ptr_free,polygon[i].caps
polygon=0

end
