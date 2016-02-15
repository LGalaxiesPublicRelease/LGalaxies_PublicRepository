;+
; NAME:
;   construct_vertex
; PURPOSE:
;   Create the structure for a vertex. 
; CALLING SEQUENCE:
;   vertex=construct_vertex([maxnvertices= ])
; INPUTS:
; OPTIONAL INPUTS:
;   maxnvertices - the maximum number of vertices allowed for any vertex
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
;   Number of caps limited to 32
; PROCEDURES CALLED:
; REVISION HISTORY:
;   01-Oct-2002  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
function construct_vertex, maxnvertices=maxnvertices

if(NOT keyword_set(maxnvertices)) then maxnvertices=15L

; make cap structure
cap1=construct_cap()

vertex={vertstr, $
        nvertices:0L, $
        weight:0.D, $
        ra:dblarr(maxnvertices), $
        dec:dblarr(maxnvertices)}

return,vertex

end
