;+
; NAME:
;   vmid
; PURPOSE:
;   Return vector within a given polygon
; CALLING SEQUENCE:
;   vec=vmid(polygon)
; INPUTS:
;   polygon - spherical polygon specification
; OPTIONAL INPUTS:
;   /allcaps - don't check whether use_caps is set, just assume all
;              caps are used
; OUTPUTS:
;   vec - [3] location of vector inside
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   16-Sep-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function vmid, polygons, allcaps=allcaps

if(n_elements(polygons) eq 0) then begin
    print, 'Syntax - vec= vmid(polygons [, /allcaps])'
    return,0
endif

; Call software
soname = filepath('libidlmangle.'+idlutils_so_ext(), $
                  root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')

vec=dblarr(3,n_elements(polygons))
for i=0L, n_elements(polygons)-1L do begin
    if(keyword_set(allcaps)) then begin
        nused_caps=polygons[i].ncaps
        used_caps=lindgen(nused_caps)
    endif else begin
        started=0
        for j=0L, polygons[i].ncaps-1L do begin
            if(is_cap_used(polygons[i].use_caps,j)) then begin
                if(started eq 0) then begin
                    used_caps=[j] 
                    started=1
                endif else begin
                    used_caps=[used_caps,j]
                endelse
            endif
        endfor
        nused_caps=n_elements(used_caps)
        if(nused_caps eq 0) then return,0
    endelse
    
    tmp_vec=dblarr(3)
    x=reform((*polygons[i].caps)[used_caps].x[*],3,nused_caps)
    cm=reform([(*polygons[i].caps)[used_caps].cm],nused_caps)
    retval = call_external(soname, 'idl_vmid', $
                           double(x), double(cm),long(nused_caps), $
                           double(tmp_vec))
    vec[*,i]=tmp_vec
endfor

return,vec
end
;------------------------------------------------------------------------------
