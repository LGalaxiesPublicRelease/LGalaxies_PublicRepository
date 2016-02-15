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
;	WAVEALIGN
;
; PURPOSE:
;	This function returns an array of arrays that has been aligned
;	to the wavelength sample of its first array.
;
; CATEGORY:
;	Pato's Miscellaneous utilities
;
; CALLING SEQUENCE:
;
;	align = WAVEALIGN(data)
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
;		aligned = WAVEALIGN(data, channel=5)
;
; MODIFICATION HISTORY:
; 	Written by:	Pato Rojo, Cornell.  2006 Jan
;			pato@astro.cornell.edu
;-


function wavealign, data, channel=channel, wavchan=wavchan, $
                    position=position, desc=desc, goodspec=goodspec, $
                    quiet=quiet, multiwav=multiwav, refcyc=refcyc

if ~ n_elements (position) then position = 0
if ~ n_elements (wavchan)  then wavchan  = 1
if ~ n_elements (channel)  then channel  = 0
if ~ n_elements (quiet)    then quiet    = 0
if ~ n_elements (multiwav) then multiwav = 0
if ~ n_elements (refcyc)   then refcyc   = 0

on_error,2

if ~ (size(data, /n_dim) eq 3 or size(data, /n_dim) eq 4) then $
  message, "Data needs to have 3 or 4 dimensions: pixel, " + $
  "cycles[, positions], channels. Instead of " + $
  string(size(data,/n_dim),format='(i0)')

nx = (size(data, /dim))[0]
nsp = (size(data, /dim))[1]
npos = 0
if size(data, /n_dim) eq 4 then begin
    npos = (size(data, /dim))[2]
    nchn = (size(data, /dim))[3]
    if ~ keyword_set(goodspec) then goodspec = bytarr(nsp, npos) + 1
endif else begin
    nchn = (size(data, /dim))[2]
    if ~ keyword_set(goodspec) then goodspec = bytarr(nsp) + 1
endelse

ret = dblarr(nx, nsp, multiwav?nsp:1)

if channel gt nchn-1 then message, $
  string("Specified channel (", channel, ") is outide boundaries: [0, ", $
         nchn-1, "]", format='(a,i0,a,i0,a)')
if wavchan gt nchn-1 then message, $
  string("Specified wavelength channel (", wavchan, $
         ") is outide boundaries: [0, ", nchn-1, "]", format='(a,i0,a,i0,a)')

if ~ quiet then $
  print, ' Interpolating ' + (keyword_set(desc)?desc:'array') + $
  ' to common wavelength... ', format='(a,$)'

if npos gt 1 then begin
    if position gt npos-1 then message, $
      "Specified position (", position, ") is outside boundaries: [0, ", $
      npos-1, "]", format='(a,i0,a,i0,a)'
    wav = data[*, *, position, wavchan]
    dat = data[*, *, position, channel]
    bs  = ~ goodspec[*, position]
endif else begin
    wav = data[*, *, wavchan]
    dat = data[*, *, channel]
    bs  = ~ goodspec
endelse

for i=0, nsp-1 do begin
    if bs[i] then continue
    c = spl_init(wav[*, i], dat[*, i], /double)

    if ~ multiwav then $
      ret[*, i] = spl_interp(wav[*,i], dat[*,i],          $
                             c, wav[*,refcyc], /double)   $
    else for ii=0, nsp-1 do begin
        if bs[ii] then continue
        ret[*, i, ii] = spl_interp(wav[*,i], dat[*,i],    $
                                   c, wav[*,ii], /double) 
    endfor
endfor

if ~ quiet then print, 'done'

return, ret

end
