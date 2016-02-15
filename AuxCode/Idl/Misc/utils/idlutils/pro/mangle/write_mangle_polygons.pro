;+
; NAME:
;   write_mangle_polygons
; PURPOSE:
;   Create a "polygon" format ascii file that mangle will understand
; CALLING SEQUENCE:
;   write_mangle_polygons, outfile, polygons [, id, weight, str, unit=]
; INPUTS:
;   outfile - output file name
;   polygons - arrays of structures (eg those made by construct_field_polygon) 
; OPTIONAL INPUTS:
;   id - array of id's for polygons (should be unique)
;   weight - arrays of weights for each polygon
;   str - area of each polygon?
;   unit - if present, use this unit instead of opening another
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
;   The format is lossy --- it only outputs "used" caps, and it throws 
;   away auxiliary information about each polygon.
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   07-Nov-2002  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
pro write_mangle_polygons, outfile, polygons, id, unit=unit

if(n_params() lt 2 or n_params() gt 3) then begin
    print,'Syntax - write_mangle_polygons, outfile, polygons [, id]'
    return
endif

if(n_elements(id) eq 0) then id=lindgen(n_elements(polygons))

if(NOT keyword_set(unit)) then $
  openw,unit,outfile,/get_lun
printf,unit,format='(%"%d polygons")',n_elements(polygons)
for i=0L, n_elements(polygons)-1L do begin
    nused_caps=0
    for j=0L, polygons[i].ncaps-1L do $
      if(is_cap_used(polygons[i].use_caps,j)) then $
         nused_caps=nused_caps+1
    printf,unit, $
      format='(%"polygon %22d ( %d caps, %20.7f weight, %20.16f str):")', $
      id[i],nused_caps,polygons[i].weight,polygons[i].str
    for j=0L, polygons[i].ncaps-1L do begin
        if(is_cap_used(polygons[i].use_caps,j)) then $
           printf,unit, $
           format='(%"%20.16f %20.16f %20.16f %20.16f")', $
           (*polygons[i].caps)[j].x[0], $
           (*polygons[i].caps)[j].x[1], $
           (*polygons[i].caps)[j].x[2], $
           (*polygons[i].caps)[j].cm
    endfor
endfor
if(NOT arg_present(unit)) then $
  free_lun,unit

end
