;+
; NAME:
;   garea
; PURPOSE:
;   Calculate area of spherical polygon by calling mangle utility garea
; CALLING SEQUENCE:
;   area=garea(polygon [, tol=, /verbose])
; INPUTS:
;   polygon - spherical polygon specification
; OPTIONAL INPUTS:
;   tol - tolerance (arcsec)
;   verbose - don't suppress garea output
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   07-Nov-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function garea, polygon, tol=tol, verbose=verbose

if(NOT keyword_set(verbose)) then $
  verbose=0L $
else $
  verbose=1L
if(n_elements(tol) eq 0) then tol=0.

for j=0L, polygon.ncaps-1L do begin
    if(is_cap_used(polygon.use_caps,j)) then begin
        if(n_elements(used_caps) eq 0) then $
          used_caps=[j] $
        else $
          used_caps=[used_caps,j]
    endif
endfor
nused_caps=n_elements(used_caps)
if(nused_caps eq 0) then return,0.D

; Call grouping software
soname = filepath('libidlmangle.'+idlutils_so_ext(), $
                  root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')
area=0.D
x=reform((*polygon.caps)[used_caps].x[*],3,nused_caps)
cm=reform([(*polygon.caps)[used_caps].cm],nused_caps)
retval = call_external(soname, 'idl_garea', $
                       double(x), double(cm), $
                       long(nused_caps), $
                       double(tol), long(verbose), area)

if(retval) then begin
    splog,'WARNING: garea.c returned error message - '+string(retval)
endif

return, area
end
;------------------------------------------------------------------------------
