pro tab_modcol,tcb,column,units,format,newname
;+
; NAME:
;	TAB_MODCOL
; PURPOSE:
;	Modify column description in a STSDAS table
;
; CALLING SEQUENCE:
;	tab_modcol, tcb, column, units, format, newname
;
; INPUTS:
;	tcb - table control block
;	column - column name or number to be modified
;
; OPTIONAL INPUTS:
;	units - string giving physical units for the column.
;		If not supplied or set to the null string
;		the units are not changed.
;	format - print format (either fortran or SPP format)
;		An spp format should be preceeded by a '%'.
;		If not supplied or set to a null string, the
;		print format for the column is not changed.
;	newname - new name for the column.  If not supplied
;		or set to a null string, the name is not
;		changed
; EXAMPLES:
;
;	change the wavelength column to WAVE with a new format
;	of 'F10.3' and columns units of ANGSTROMS.
;
;	   tab_modcol,tcb,'wavelength','ANGSTROMS','F10.3','WAVE'
;
;	Change to print format of column 3 to spp format
;	20.10e
;	   tab_modcol,tcb,3,'','%20.10e'
; HISTORY:
;	version 1  D. Lindler   Apr 1989
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-----------------------------------------------------------------------------
;
; set defaults for missing parameters
;
	if n_params(0) lt 3 then units=''
	if n_params(0) lt 4 then format=''
	if n_params(0) lt 5 then newname=''
;
; check that the column exists and get its number
;
	tab_col,tcb,column
	colnum=!err
	if colnum lt 1 then begin
		print,'TAB_MODCOL-- column '+strtrim(column)+' does not exist'
		retall
	endif
;
; update name
;
	col_desc=byte(tcb[*,colnum],0,64)
	if strtrim(newname) ne '' then begin
		name=bytarr(19)
		name[0]=byte(strtrim(newname))
		col_desc[16]=name
	endif
;
; update column units
;
	if units ne '' then begin
		u=bytarr(19)
		u[0]=byte(units)
		col_desc[36]=u
	endif
;
; update print format
;
	if strtrim(format) ne '' then begin
		if strmid(format,0,1) eq '%' then $	;spp format?
			sppformat=strmid(format,1,strlen(format)-1) $
			else tab_fortospp,format,sppformat
		f = bytarr(8)
		f[0]=byte(sppformat)
		col_desc[56]=f
	endif
;
; place new column description into the table control block
;
	tcb[0,colnum]=long(col_desc,0,16)
return
end	
