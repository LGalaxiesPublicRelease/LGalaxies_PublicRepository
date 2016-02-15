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
;	WAVEDEALIGN
;
; PURPOSE:
;	This function returns an array of arrays with a common sample
;	to one of the given wavelength sample
;
; CATEGORY:
;	Pato's data reduction potpourri.
;
; CALLING SEQUENCE:
;
;	align = WAVEDEALIGN(data)
;
; INPUTS:
;	Parm1:	Initial caps
;
; OPTIONAL INPUTS:
;	Parm2:	Initial caps
;	
; KEYWORD PARAMETERS:
;	MULTIWAV: It signals that the returned array should have three
;                 dimensions where the third indicates the wavelength 
;	          basis.
;
; OUTPUTS:
;	This function returns an array with the appropiate
;	cycles[,positions] aligned according to the firstr solution.
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
;		aligned = WAVEDEALIGN(data, channel=5)
;
; MODIFICATION HISTORY:
; 	Written by:	Pato Rojo, Cornell.  2006 Jan
;			pato@astro.cornell.edu
;-

function wavedealign, wav, data,                    $
                      refpos=refpos, refwav=refwav, $
                      desc=desc, goodspec=goodspec, $
                      quiet=quiet

if ~ keyword_set(refwav) then $
  refwav = wav[*,n_elements(position)?position:0]

if ~ n_elements (quiet)    then quiet    = 0
if ~ n_elements (multiwav) then multiwav = 0
if ~ n_elements (refcyc)   then refcyc   = 0

on_error,2

if ~ (size(data, /n_dim) eq 2) then $
  message, "Data needs to have 3 dimensions: pixel, " + $
  "cycles. Instead of " + $
  string(size(data,/n_dim),format='(i0)')

nx = (size(data, /dim))[0]
nsp = (size(data, /dim))[1]
if ~ keyword_set(goodspec) then goodspec = bytarr(nsp) + 1

ret = dblarr(nx, nsp)

if ~ quiet then $
  print, ' Interpolating ' + (keyword_set(desc)?desc:'array') + $
  ' to each wavelength... ', format='(a,$)'

for i=0, nsp-1 do begin
    if ~ goodspec[i] then continue
    c = spl_init(refwav, data[*, i], /double)
    ret[*, i] = spl_interp(refwav, data[*,i], c, wav[*,i], /double)
endfor

if ~ quiet then print, 'done'

return, ret

end
