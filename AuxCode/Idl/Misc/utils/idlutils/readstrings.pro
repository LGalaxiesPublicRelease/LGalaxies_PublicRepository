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
;	READSTRINGS
;
; PURPOSE:
;	This function  reads lines from a file
;
; CATEGORY:
;	Pato's data miscellaneous utils
;
; CALLING SEQUENCE:
;
;	strings = readstrings(Filename)
;
; INPUTS:
;	Parm1:	Initial caps
;
; KEYWORD PARAMETERS:
;	KEY1:	ALL CAPS!
;
; PROCEDURE:
;	Algorithm description.
;
; EXAMPLE:
;	Describe example here
;
;		F = READSTRINGS('data.txt')
;
; MODIFICATION HISTORY:
; 	Written by:	Pato Rojo, Cornell.  2006
;			pato@astro.cornell.edu
;-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function readstrings, file, nlines=nlines

; read an ASCII file into an array of strings, one string per line

; file		filename
; nlines	(returned) number of lines in file

line   = 'this is a string'
nlines = 0L

; first count the lines...hope this isn't too slow
openr, fnum, file, /get_lun
repeat begin
  readf, fnum, line
  nlines = nlines + 1
end until eof(fnum)
free_lun, fnum

asctab = strarr(nlines)

openr, fnum, file, /get_lun
for i = 0L, nlines-1, 1 do begin
  readf, fnum, line
  asctab[i] = line
endfor
free_lun, fnum

return, asctab
end

