pro table_sort,name,column,name_out
;+
; NAME:
;	TABLE_SORT
; PURPOSE:
;	Procedure to sort an STSDAS table by the specified column
;
; CALLING SEQUENCE:
;	table_sort, name, column, [ name_out ]
;
; INPUTS:
;	name - table name
;	column - column to sort on
;
; OPTIONAL INPUTS:
;	name_out - output table name.  If not supplied, input name
;		is used.
;
; HISTORY:
;	version 1  D. Lindler  MAY 1989
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-----------------------------------------------------------------------
	if n_params(0) lt 3 then out_name=name
	tab_read,name,tcb,tab,h
	tab_sort,column,tcb,tab
	tab_write,out_name,tcb,tab,h
	end
