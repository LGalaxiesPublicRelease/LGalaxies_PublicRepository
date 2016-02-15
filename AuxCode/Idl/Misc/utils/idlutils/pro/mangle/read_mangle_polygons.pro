;+
; NAME:
;   read_mangle_polygons
; PURPOSE:
;   Read a "polygon" format ascii file written by mangle, and return
;   in the IDL structure format
; CALLING SEQUENCE:
;   read_mangle_polygons, infile, polygons, id [, unit=]
; INPUTS:
;   infile - input file name
; OPTIONAL INPUTS:
;   unit - if present, read from given unit
; OUTPUTS:
;   polygons - arrays of structures (eg those made by construct_field_polygon) 
;   id - array of id's for polygons (should be unique)
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   30-Nov-2002  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
pro read_mangle_polygons, infile, polygons, id, unit=unit, $
                          allow_doubles=allow_doubles

if(n_params() lt 2) then begin
    print,'Syntax - read_mangle_polygons, infile, polygons [,id, unit=]'
    return
endif

if(NOT keyword_set(unit)) then $
  openr,unit,infile,/get_lun
npoly=0L
readf,unit, format='(i,"polygons")',npoly
tmp_line=''
id=lon64arr(npoly)
polygons=replicate(construct_polygon(),npoly)
for i=0L, npoly-1L do begin
    readf,unit, tmp_line
    tmp_words=strsplit(tmp_line,/extract)
    id[i]=long64(tmp_words[1])
    ptr_free,polygons[i].caps
    polygons[i].weight=double(tmp_words[5])
    polygons[i].str=double(tmp_words[7])
    polygons[i].ncaps=long(tmp_words[3])
    polygons[i].caps=ptr_new(replicate(construct_cap(),polygons[i].ncaps))
    for j=0L, polygons[i].ncaps-1L do begin
       readf,unit,tmp_line
       tmp_words=strsplit(tmp_line,/extract)
       (*(polygons[i].caps))[j].x[0:2]=double(tmp_words[0:2])
       (*(polygons[i].caps))[j].cm=double(tmp_words[3])
    endfor
    set_use_caps,polygons[i],lindgen(polygons[i].ncaps), use_caps=use_caps, $
      allow_doubles=allow_doubles
    polygons[i].use_caps=use_caps
endfor
if(NOT arg_present(unit)) then $
  free_lun,unit

end
