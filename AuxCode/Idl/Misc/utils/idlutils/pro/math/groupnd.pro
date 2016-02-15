;+
; NAME:
;   groupnd
; PURPOSE:
;   group a set of points in N dimensions
; CALLING SEQUENCE:
;   ingroup= groupnd(x, distance [, nextgroup=, multgroup=, $
;                    firstgroup= ])
; INPUTS:
;   x - [M,N] positions in M-dimensions
;   distance - link distance
; OUTPUTS:
;   ingroup    - group number of each object (N-dimensional array);
;                -1 if no groups
;   multgroup  - multiplicity of each group 
;   firstgroup - first member of each group 
;   nextgroup  - index of next member of group for each object
; REVISION HISTORY:
;   28-Nov-2006  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function groupnd, x, distance, nextgroup=nextgroup, multgroup=multgroup, $
                  firstgroup=firstgroup, nd=nd

if((size(x))[0] eq 1) then begin
    if(NOT keyword_set(nd)) then begin
        mm=1
        nn=n_elements(x)
        x=reform(x, 1, nn)
    endif else begin
        mm=nd
        nn=n_elements(x)/mm
        x=reform(x, mm, nn)
    endelse 
endif else begin
    mm=(size(x,/dim))[0]
    nn=(size(x,/dim))[1]
endelse

matchnd, x, x, distance, m1=m1, m2=m2, d12=d12, nmatch=nmatch, nd=nd, $
  maxmatch=0

matches=ptrarr(nn)
isort=sort(m1)
iuniq=uniq(m1[isort])
istart=0L
for i=0L, n_elements(iuniq)-1L do begin
    iend=iuniq[i]
    icurr=isort[istart:iend]
    matches[m1[icurr[0]]]=ptr_new(m2[icurr])
    istart=iend+1L
endfor

group_on_matches, matches, first=firstgroup, next=nextgroup, $
  mult=multgroup, in=ingroup

return, ingroup

end
