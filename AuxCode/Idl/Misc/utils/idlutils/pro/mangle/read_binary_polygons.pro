;+
; NAME:
;   read_binary_polygons
; PURPOSE:
;   Read a "polygon" format binary file written by mangle, and return
;   in the IDL structure format
; CALLING SEQUENCE:
;   read_binary_polygons, infile, polygons, id [, unit=]
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
;   27-Sep-2003  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
pro read_binary_polygons, infile, polygons, id, unit=unit, $
                          allow_doubles=allow_doubles

if(n_params() lt 2) then begin
    print,'Syntax - read_binary_polygons, infile, polygons [,id, unit=]'
    return
endif

if(NOT keyword_set(unit)) then $
  openr,unit,infile,/get_lun
npoly=0L
readu,unit, npoly
tmp_line=''
id=lon64arr(npoly)
polygons=replicate(construct_polygon(),npoly)
tmp_id=0L
tmp_weight=0.D
tmp_np=0L
tmp_str=0.D
for i=0L, npoly-1L do begin
    readu,unit, tmp_id
    readu,unit, tmp_np
    readu,unit, tmp_weight
    readu,unit, tmp_str
    id[i]=tmp_id
    ptr_free,polygons[i].caps
    polygons[i].ncaps=tmp_np
    polygons[i].weight=tmp_weight
    polygons[i].str=tmp_str
    polygons[i].caps=ptr_new(replicate(construct_cap(),polygons[i].ncaps))
    tmp_x=dblarr(3,tmp_np)
    readu,unit,tmp_x
    tmp_cm=dblarr(tmp_np)
    readu,unit,tmp_cm
    for j=0L, polygons[i].ncaps-1L do begin
        (*(polygons[i].caps))[j].x[0:2]=tmp_x[*,j]
        (*(polygons[i].caps))[j].cm=tmp_cm[j]
    endfor
    set_use_caps,polygons[i],lindgen(polygons[i].ncaps), use_caps=use_caps, $
      allow_doubles=allow_doubles
    polygons[i].use_caps=use_caps
endfor
if(NOT arg_present(unit)) then $
  free_lun,unit

end
