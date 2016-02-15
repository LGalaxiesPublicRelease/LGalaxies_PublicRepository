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
;	FITSFIELD
;
; PURPOSE:
;	This function return an array with the value of a specified
;	FITS field.
;
; CATEGORY:
;	Pato's data reduction potpourri
;
; CALLING SEQUENCE:
;
;	Result = FITSFIELD(Header, Name, Errormessage)
;
; INPUTS:
;	Header:   Header array
;       Name:     Name of the field
;
; OPTIONAL INPUTS:
;       Errormessage: If present. The code stop with that error
;                     messsage when the field is not found.
;
; OUTPUTS:
;	This function returns an array of strings with the appropriate
;	field.
;
; EXAMPLE:
;	Describe example here
;
;           time = double(fitsfield(header, 'EXPTIME', 'EXPTIME not found'))
;
; MODIFICATION HISTORY:
; 	Written by:	Pato Rojo, Cornell.  2006
;			pato@astro.cornell.edu
;-

function fitsfield, headarr, name, errormessage

if n_params() lt 3 then errormessage = 0

nfiles = size(headarr, /n_dim) eq 2? (size(headarr, /dim))[1]: 1
val = strarr(nfiles)

len = strlen(name)
hierarch = strmid(name, 0, 8) eq 'HIERARCH'

;check that name is right
if len lt 8 then fname = strupcase(name) + strjoin(replicate(" ", 8-len)) $
else if len eq 8 or hierarch then fname = strupcase(name) $
else begin
    if keyword_set(errormessage) then begin
        print, "Field name '" + name + "' is not a valid fits field's name"
        message, errormessage
    endif else return, -1
endelse

;add equal sign unles it is a comment or an ESO HIERARCH
;keyword. Mosty to avoid mismatching.
if strmid(fname, 0, 7) ne 'COMMENT' and ~ hierarch then begin
    case strlen(stregex(name, "= ?$", /extr)) of
        0: fname = fname + "= "
        1: fname = fname + " "
    endcase
endif

len = strlen(fname)

;look for the value
for i=0, nfiles-1 do begin
    idx = where(strpos(headarr[*, i], fname) eq 0)
    if n_elements(idx) gt 1 or idx[0] eq -1 then begin
        if keyword_set(errormessage) then begin
            if idx[0] eq -1 then print, "Field name '" + name + $
              "' was not found" $
            else print, "There is more than one match for field name '" + $
              name + "'. Aborting"
            message, errormessage
        endif else begin
            val[i] = -1
            continue
        endelse
    endif
;remove equal sign if it is an hierarch keyword
    fc = strmid(headarr[idx[0], i], len)
    if hierarch then fc = strmid(fc, strlen(stregex(fc, "^ *=", /extr)))
    fc = strtrim(fc, 1)
    mc = strpos(fc, "/")
    cand = strmid(fc, 0, mc eq -1? 80: mc)
;    cand = stregex(fc, "^ *([^ ]+) *(/.*)?$", /extr, /subs)
    if strlen(cand) gt 0 then val[i] = cand else val[i] = -1
endfor

return, val

end
