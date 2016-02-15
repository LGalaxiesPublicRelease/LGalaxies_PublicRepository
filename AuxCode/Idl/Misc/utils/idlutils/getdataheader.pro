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
;	GETDATAHEADER
;
; PURPOSE:
;	This function reads data and particular header fields from
;	FITS files.
;
; CATEGORY:
;	Pato's data reduction potpourri
;
; CALLING SEQUENCE:
;
;         GETDATAHEADER, Files, Data, Head, Type, Res, Nx, Ny
;
; INPUTS:
;	Files:  FITS files to read
;       Head:   Header that want to be read
;       Type:   Type of the header to be read. 0 for double, 1 for
;               boolean.
;
; KEYWORD PARAMETERS:
;	KEY1:	ALL CAPS!
;
; OUTPUTS:
;       Data:   Returns data array
;       Res:    Returns header values, is always a double but it can
;               have boolean meaning
;       Nx:     Returns size in X
;       Ny:     Returns size in Y
;
; PROCEDURE:
;	Algorithm description.
;
; EXAMPLE:
;	Describe example here
;
;
; MODIFICATION HISTORY:
; 	Written by:	Pato Rojo, Cornell.  2006
;			pato@astro.cornell.edu
;-

pro getdataheader, files, data, head, type, res, nx, ny, $
                   header=oh, dot=dot, instrument=instrument


sd = keyword_set(dot)

nh = n_elements(head);
nf = n_elements(files);
res = dblarr(nh, nf);

data = readfits(files[0], /silent)
nx = (size(data, /dim))[0]
ny = (size(data, /dim))[1]

data = fltarr(nx, ny, nf)
oh = strarr(nf, 600)

if ~ keyword_set(instrument) then instrument = 0

for i = 0, nf-1 do begin
    if sd then print, ".", format='(a,$)'
    data[*,*,i] = readfits(preprocfits(files[i], instrument, silent=sd), $
                           h, /silent)
    oh[i,0:n_elements(h)-1] = h

    for j = 0, nh-1 do begin
        mjd = where(strpos(h, head[j]) ne -1, njd)
        if (njd le 0) then                                    $
          message, "File " + files[i] + " doesn't have a '" + $
          head[j] + "' header"
        case type[j] of
            0: begin
                sjd = stregex(h[mjd], '[a-z]+ *= *([0-9]+(\.[0-9]+)?)',  $
                              /fold_case, /subexpr, /extract)
                res[j,i] = double(sjd[1])
            end
            1: begin
                sjd = stregex(h[mjd], '[a-z]*= *(T|F)',             $
                              /fold_case, /subexpr, /extract)
                res[j,i] = double(sjd[1] eq 'T')
            end
        endcase
              
    endfor

endfor
if sd then print, " "

dummy = preprocfits('whatever', 'isaac', /delete)

end
