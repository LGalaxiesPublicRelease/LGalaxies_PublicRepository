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
;	ARRAYFROMFILE
;
; PURPOSE:
;	This function returns a double precision array with the contents
;	of column separated data file
;
; CATEGORY:
;	Pato's miscellaneous utilities.
;
; CALLING SEQUENCE:
;
;	data = ARRAYFROMFILE(File)
;
; INPUTS:
;	Parm1:	Initial caps
;
; OPTIONAL INPUTS:
;	Parm2:	Initial caps
;	
; KEYWORD PARAMETERS:
;	KEY1:	ALL CAPS!
;
; OUTPUTS:
;	This function returns 
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;	BLOCK1:	All caps
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;	Algorithm description.
;
; EXAMPLE:
;	Describe example here
;
;		F = PICKFILE(/READ, FILTER = ['pro', 'dat'])
;
; MODIFICATION HISTORY:
; 	Written by:	Pato Rojo, Cornell.  2006
;			pato@astro.cornell.edu
;       Based on code by Joseph Harrington.
;-


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function decomment_af, asctab, comchar=comchar, nlines=nlines, maxtok=maxtok

; remove comments, leading/trailing blanks from lines in a string array,
; and remove blank lines

; asctab	an array of strings ("ASCII table")
; comchar	comment character, default "#"
; nlines	(returned) number of lines left in array
; maxtok	(returned) maximum number of tokens on a line
;		Not computed unless required.  MUST BE DEFINED BEFORE CALL!

if not keyword_set(comchar) then begin
  comchar = '#'
endif

nl    = n_elements(asctab)
out   = strarr(nl)
kdmt  = n_elements(maxtok) gt 0
nlines = 0L

for i = 0L, nl-1, 1 do begin
  line = asctab[i]
  n = strpos(line, comchar)
  if n ne -1 then begin
    line = strmid(line, 0, n)
  endif
  line = strtrim(line, 2)
  if strlen(line) gt 0 then begin
    out[nlines] = line
    nlines = nlines + 1
    if kdmt then begin
      tok = strsplit(line)
      maxtok = maxtok > n_elements(tok)
    endif
  endif
endfor

return, out[0:nlines-1]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function parsetab_af, table, comchar=comchar

; read an ascii table into a 2D array of strings, one per word

maxtok = 0
at   = decomment_af(table, nlines=nlines, maxtok=maxtok, comchar=comchar)
ptab = strarr(maxtok, nlines)

for j = 0L, nlines-1, 1 do begin
  tok = strsplit(at[j], /extract)
  ptab[0,j] = tok[*]
endfor

return, transpose(ptab)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function arrayfromfile, file, nlines=nlines, comchar=comchar

; read an ASCII file into an array of doubles

; file		filename
; nlines	(returned) number of lines in file

return, double(parsetab_af(readstrings(file, nlines=nlines), $
                           comchar=comchar))

end
