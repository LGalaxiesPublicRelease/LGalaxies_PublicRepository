;+
; NAME:
;   read_fits_polygons
; PURPOSE:
;   Read a "polygon" format fits file 
; CALLING SEQUENCE:
;   read_fits_polygons, infile, polygons
; INPUTS:
;   infile - input file name
; OPTIONAL INPUTS:
; OUTPUTS:
;   polygons - arrays of structures (eg those made by construct_field_polygon) 
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
;   The main point of this is to extract the xcaps and cmcaps columns
;   and replace them with caps.x and caps.cm
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   30-Nov-2002  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
pro read_fits_polygons, infile, polygons, _EXTRA=extra_for_mrdfits

if(n_params() ne 2) then begin
    print, 'Syntax - read_fits_polygons, infile, polygons'
    return
endif

inpoly=mrdfits(infile,1,/unsigned,/silent,_EXTRA=extra_for_mrdfits)
intags=tag_names(inpoly)
for i=0L, n_elements(intags)-1L do begin
    if(intags[i] ne 'XCAPS' and $
       intags[i] ne 'CMCAPS') then begin
        if(n_tags(polygon1) gt 0) then $
          polygon1=create_struct(polygon1,intags[i], inpoly[0].(i)) $
        else $
          polygon1=create_struct(intags[i], inpoly[0].(i)) 
    endif 
endfor
cap1=construct_cap()
polygon1=create_struct(polygon1,'caps',ptr_new(0))

polygons=replicate(polygon1,n_elements(inpoly))
struct_assign,inpoly,polygons,/nozero

for i=0L, n_elements(inpoly)-1L do begin
    polygons[i].caps=ptr_new(replicate(construct_cap(),inpoly[i].ncaps))
    (*polygons[i].caps)[0:inpoly[i].ncaps-1].x= $
      inpoly[i].xcaps[*,0:inpoly[i].ncaps-1]
    (*polygons[i].caps)[0:inpoly[i].ncaps-1].cm= $
      inpoly[i].cmcaps[0:inpoly[i].ncaps-1]
endfor

end
