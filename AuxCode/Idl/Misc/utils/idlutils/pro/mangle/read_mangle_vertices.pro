;+
; NAME:
;   read_mangle_vertices
; PURPOSE:
;   Read in a set of vertices, mangle style
; CALLING SEQUENCE:
;   read_mangle_vertices, infile, vertices, id
; INPUTS:
;   infile - input file name
; OPTIONAL INPUTS:
; OUTPUTS:
;   vertices - arrays of structures (eg those made by
;              construct_field_vertex) 
;   id - array of id's for vertices (should be unique)
;   weight - arrays of weights for each vertex
;   str - area of each vertex?
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   30-Nov-2002  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
pro read_mangle_vertices, infile, vertices, id, maxnvertices=maxnvertices

if(n_params() lt 2) then begin
print, 'Syntax - read_mangle_vertices, infile, vertices [, id]'
return
endif

if(NOT keyword_set(maxnvertices)) then maxnvertices=15

openr,unit,infile,/get_lun
nvertex=0L
readf,unit, format='(i,"polygons")',nvertex
tmp_line=''
readf,unit, tmp_line
id=lon64arr(nvertex)
vertices=replicate(construct_vertex(maxnvertices=maxnvertices),nvertex)
indx0=2*lindgen(maxnvertices)
indx1=indx0+1L
for i=0L, nvertex-1L do begin
    readf,unit, tmp_line
    tmp_words=strsplit(tmp_line,/extract)
    id[i]=long64(tmp_words[1])
    tmp_vertex=construct_vertex(maxnvertices=maxnvertices)
    tmp_vertex.nvertices=long(tmp_words[3])
    tmp_vertex.weight=double(tmp_words[5])
    if(tmp_vertex.nvertices gt 0) then begin
        readf,unit,tmp_line
        tmp_words=strsplit(tmp_line,/extract)
        tmp_vertex.ra=double(tmp_words[indx0[0L:tmp_vertex.nvertices-1L]])
        tmp_vertex.dec=double(tmp_words[indx1[0L:tmp_vertex.nvertices-1L]])
    endif
    vertices[i]=tmp_vertex
endfor
free_lun,unit

end
