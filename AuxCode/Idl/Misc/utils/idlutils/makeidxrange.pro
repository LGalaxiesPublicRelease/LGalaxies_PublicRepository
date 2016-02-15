;;
; 
; Copyright (C) 2006 Patricio Rojo
; 
; This program is free software; you can redistribute it and/or
; modify it under the terms of the GNU General Public License
; as published by the Free Software Foundation; either version 2
; of the License, or (at your option) any later version.
; 
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
; 
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
;   MA  02110-1301, USA.
; 
;;

;+
; NAME:
;	MAKEIDXRANGE
;
; PURPOSE:
;	This function makes an index range from a value range
;
; CATEGORY:
;	Pato's miscellaneous utilities.
;
; CALLING SEQUENCE:
;
;	MAKEIDXRANGE, Lims, Xaxis
;
; INPUTS:
;	Lims:        Limit of the range
;       Xaxis:       X-axis values.
;	
; KEYWORD PARAMETERS:
;	MAPDEG:      Degree of X-axis mapping.
;
; OUTPUTS:
;	This function returns 
;
; PROCEDURE:
;	Algorithm description.
;
; EXAMPLE:
;	Describe example here
;
;		MAKEIDXRANGE, Lims, Xaxis
;
; MODIFICATION HISTORY:
; 	Written by:	Pato Rojo, Cornell.  2006
;			pato@astro.cornell.edu
;-

function polyn_locidxrng, c, x

ret = 0
for i=n_elements(c)-1, 0, -1 do ret = x*ret + c[i]

return, ret

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function makeidxrange, lims, xaxis, mapdeg=mapdeg

if not keyword_set(mapdeg) then  mapdeg = 2
if n_params() eq 1 then idxlim = long(lims+0.499) $
else begin
    coef = poly_fit(xaxis, indgen(n_elements(xaxis)), mapdeg)
    idxlim = long(polyn_locidxrng(coef, lims) + 0.499)
    if size(idxlim, /n_dim) eq 1 then idxlim = reform(idxlim, 2, 1)
endelse

if (size(idxlim, /n_dim)) ne 2 then $
  message, "Range boundary is not well set in makeidxrange"
n = (size(idxlim, /dim))[1]
for i=0, n-1 do begin
    mx = max(idxlim[*,i], min=mn)
    idxlim[*, i] = [mn, mx]
endfor

;number of indices per pair
ni = idxlim[1,*] - idxlim[0,*] + 1

;if there is more than one pair
if (size(idxlim, /dim))[1] gt 1 then begin
    tx = lindgen((size(idxlim,/dim))[1]-1)
    idx = lindgen(total(ni))
    acum = 0
    diff = [idxlim[0,0], reform(idxlim[0,tx+1] - idxlim[1,tx]) - 1]

    for i=0, n_elements(tx) do begin
        idx[acum:*] += diff[i]
        acum += ni[i]
    endfor
;if there is only one pair
endif else idx = lindgen(ni[0]) + idxlim[0, 0]

return, idx[uniq(idx[where(idx ge 0)], sort(idx))]

end
