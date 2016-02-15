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
; 
;;

;+
; NAME:
;	MEDIANGOODP
;
; PURPOSE:
;       Make a good pixel mask using median method
;
; CATEGORY:
;
; CALLING SEQUENCE:
;
;	Result = mediangoodp(data, sig, rad, [/showk])
;
; INPUTS:
;       Data:   Array of two dimensional data as a function of time
;       Sigma:  Number of sigma for rejection of bad pixels
;       Rad:    Radius of area in which to process bad pixels
;
; KEYWORD PARAMETERS:
;       showk:  Whether to plot the difference and pixel mask
;
; OUTPUTS:
;	This function returns a good pixel mask (1 for good, 0 for
;	bad)
;
; SIDE EFFECTS:
;        Occasional nausea
;
; RESTRICTIONS:
;
; PROCEDURE:
;	Checks for bad pixels and store them in a 3D integer
;	array, bad pixels are found by looking for values that go 'sig'
;	times over the median in a box of radius 'rad'.
;
; EXAMPLE:
;
;		goodp = mediangoodp(data,5,15)
;
; MODIFICATION HISTORY:
; 	Written by:	Pato Rojo, Cornell.  2004-2005
;			pato@astro.cornell.edu
;-
function mediangoodp, data, sig, rad, showk=showk, $
                     databorders=databorders

nx = (size(data, /dim))[0]
ny = (size(data, /dim))[1]
nf = (size(data, /n_dim) eq 3)?(size(data, /dim))[2]:1

if keyword_set(databorders) then begin
    x0 = databorders[0]
    x1 = nx-databorders[2]-1
    y0 = databorders[1]
    y1 = ny-databorders[3]-1
endif else begin
    x0 = 0
    x1 = nx-1
    y0 = 0
    y1 = ny-1
endelse

show = keyword_set(showk)

print, "Computing good pixel mask", format='(a,$)'
goodp = bytarr(nx, ny, nf)
for i = 0, nf-1 do begin
    med = median(reform(data[*,*,i]), rad, /double)
    diff = reform(data[*,*,i]) - med

    sd = stddev(diff[x0:x1,y0:y1])
    goodp[*,*,i] = abs(diff) lt sig*sd

    if show then begin
        erase
        !p.multi = [2,2,1]
        shade_surf, diff
        shade_surf, goodp[*,*,i]
        !p.multi = 0
    endif
    print, ".", format='(a,$)'
endfor
print, " "

return, goodp
            
end
