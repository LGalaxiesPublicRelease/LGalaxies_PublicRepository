pro tab_del,tcb,tab,rows
;+
; NAME:
;	TAB_DEL 
; PURPOSE:
;	Delete specified row(s) from an STSDAS table
;
; CALLING SEQUENCE:
;	tab_del, tcb, tab, rows
;
; INPUT/OUTPUTS
;	tcb - table control block
;	tab - table array
;
; OPTIONAL INPUTS:
;	rows - row (scalar) or rows(vector) to delete from the table
;		If not supplied all rows are deleted.
;
; HISTORY:
;	version 1  D. Lindler  April 1989
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-----------------------------------------------------------------------
;
; If rows not supplied, delete all rows
;
	if n_params(0) lt 3 then begin
		tcb[3]=0
		return
	endif
;
; make rows into a vector
;
	ndel=n_elements(rows)			;number of rows to delete
	r=lonarr(ndel)+rows
	nrows=tcb[3]				;number of rows in the table
	if (max(r) ge nrows) or (min(r) lt 0) then begin
		print,'TAB_DEL-- Invalid row number specified'
		retall
	endif
;
; create a mask of rows to keep
;
	mask=replicate(1,nrows)
	mask[r]=0			;flag lines to delete
	keep=where(mask)		;vector of rows to keep
	nkeep=!err			;number of rows to keep
	tcb[3]=nkeep>0
	if nkeep lt 1 then return	;all rows deleted?
;
; compress rows kept
;
	pos=0				;output position
	for i=0,nkeep-1 do begin
		tab[0,pos]=tab[*,keep[i]]
		pos=pos+1
	endfor
return
end
