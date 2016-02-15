pro table_delete,name,rows,outname
;+
; NAME:
;	TABLE_DELETE 
; PURPOSE:
;	Delete specified rows from an STSDAS table
;
; CALLING SEQUENCE:
;	table_delete, name, rows, [ outname ]
;
; INPUT:
;	name - table name
;	rows - row (scalar) or rows(vector) to delete from the table
;
; OPTIONAL OUTPUT:
;	outname - output table name, if not supplied the input name
;		is used
;
; HISTORY:
;	version 1  D. Lindler  April 1989
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-----------------------------------------------------------------------
;
	if n_params(0) lt 2 then begin
		print,'TABLE_DELETE-- at least two parameters required'
		retall
	endif
	if n_params(0) lt 3 then outname=name
	tab_read,name,tcb,tab,h		;read table
	tab_del,tcb,tab,rows		;delete rows
	tab_write,outname,tcb,tab,h	;write table
	return
	end
