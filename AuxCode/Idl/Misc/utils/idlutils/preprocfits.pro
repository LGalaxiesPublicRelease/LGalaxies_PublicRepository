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
;	PREPROCFITS
;
; PURPOSE:
;	This function should be run before processing any FITS.
;
; CATEGORY:
;	Pato's data reduction potpourri
;
; CALLING SEQUENCE:
;
;       newfilename =  PREPROCFITS(Filename, Instrument)
;
; INPUTS:
;	Filename:  FITS filename
;       Instrument: Instrument name. Right now only 'ISAAC' is available.
;
; KEYWORDS:
;
; OUTPUTS:
;	This function returns the filename of a FITS file processed by
;	this algorithm
;
; OPTIONAL OUTPUTS:
;
; PROCEDURE:
;	Algorithm description.
;
; EXAMPLE:
;	Describe example here
;
;		a = READFITS(PREPROCFITS('data.fits', 'isaac'))
;
; MODIFICATION HISTORY:
; 	Written by:	Pato Rojo, Cornell.  2006
;			pato@astro.cornell.edu
;-

function preprocfits, filename, instrument, delete=delete, silent=silent


dir = file_dirname(filename, /mark)
tmpfile = file_basename(filename)

if instrument eq 'isaac' then begin
    tmpfile = 'isaacdeghosted.tmp.fits'
    if keyword_set(delete) then begin
        if ~ keyword_set(silent) then $
          print, "Removing temporary file used for preprocessing"
        spawn, 'chmod +rw ' + tmpfile
        file_delete, tmpfile
        return, filename
    endif
    if ~ keyword_set(silent) then $
      print,"    => Preprocessing ISAAC's FITS file"
    if (file_test(tmpfile)) then $
      spawn, 'chmod +rw ' + tmpfile
    file_copy, filename, tmpfile, /overwrite
    spawn, 'chmod +rw ' + tmpfile
    spawn, "isaacp ghost " + tmpfile, stdin, stdout
endif

return, tmpfile

end
