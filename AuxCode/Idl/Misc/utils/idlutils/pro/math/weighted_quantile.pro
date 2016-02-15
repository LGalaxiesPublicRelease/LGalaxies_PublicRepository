;+
; NAME:
;   weighted_quantile
; PURPOSE:
;   given a set of values and weights, returns weighted quantile(s)
; CALLING SEQUENCE:
;   quantile= weighted_quantile(values,weights,quant=0.25)
; REVISION HISTORY
;   2002-07-ish  written - Blanton
;-
function weighted_quantile,values,weights,quant=quant

if not keyword_set(weights) then weights= dblarr(n_elements(values))+1
if(n_elements(values) le 1) then return,0.*quant+values[0]
if(n_elements(quant) eq 0) then quant=[0.5]

isort=sort(values)
svalues=values[isort]
sweights=weights[isort]
scum=total(sweights,/cumulative,/double)
scum=scum/scum[n_elements(scum)-1L]
j=lindgen(n_elements(scum)-1L)
jp1=j+1L

quantile=dblarr(n_elements(quant))
for iquant=0L, n_elements(quant)-1L do begin
    ipos=where(scum[j] le quant[iquant] and $
               scum[jp1] gt quant[iquant],ispos)
    quantile[iquant]=svalues[n_elements(svalues)-1L]
    if(scum[0] gt quant[iquant]) then quantile[iquant]=svalues[0]
    if(ispos gt 0) then quantile[iquant]=svalues[ipos[0]+1]
endfor

if n_elements(quant) EQ 1 then quantile=quantile[0]

return,quantile

end
