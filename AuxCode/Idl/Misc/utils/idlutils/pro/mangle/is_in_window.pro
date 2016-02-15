;+
; NAME:
;   is_in_window
; PURPOSE:
;   Is an xyz (or radec) position in any of a given list of polygons?
; CALLING SEQUENCE:
;   result=is_in_window(xyz=, ra=, dec= , polygons)
; INPUTS:
;   polygons - polygons with caps to check
; OPTIONAL INPUTS:
;   ra - [N] set of ra values
;   dec - [N] set of dec values
;   xyz - [3,N] xyz value(s) (overrides ra and dec)
;   ncaps - override polygon.ncaps (if ncaps < polygon.ncaps)
; OUTPUTS:
;   result - [N] 1 if in window, 0 otherwise
; OPTIONAL OUTPUTS:
;   in_polygon - [N] which polygon each ra,dec is in (-1 if none)
; COMMENTS:
;   Either ra and dec, or xyz must be set; xyz overrides ra and dec
; REVISION HISTORY:
;   01-Oct-2002  Written by MRB (NYU)
;-
;------------------------------------------------------------------------------
function is_in_window, polygons, ra=ra, dec=dec, ncaps=ncaps, $
                       in_polygon=in_polygon

nxyz=n_elements(ra) 
in_window=lonarr(nxyz)
in_polygon=lonarr(nxyz)-1L
curr_polygon=0L
while(curr_polygon lt n_elements(polygons)) do begin
    indx_not_in=where(in_polygon eq -1L,count_not_in)
    if(count_not_in gt 0) then begin
        indx_in_curr_polygon= $
          where(is_in_polygon(ra=ra[indx_not_in],dec=dec[indx_not_in], $
                              polygons[curr_polygon], ncaps=ncaps), $
                count_in_curr_polygon)
        if(count_in_curr_polygon gt 0) then $
          in_polygon[indx_not_in[indx_in_curr_polygon]]=curr_polygon
    endif
    curr_polygon=curr_polygon+1L
endwhile
in_window=(in_polygon ge 0L)
return,in_window

end
