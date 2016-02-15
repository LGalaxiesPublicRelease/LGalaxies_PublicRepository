pro table_print,name,columns,row1,row2
;
;+
; NAME:
;	TABLE_PRINT
; PURPOSE:
;	Routine to print an stsdas table.
;
; CALLING SEQUENCE:
;	table_print, name, columns, row1, row2
;
; INPUTS:
;	name - table name
; 
; OPTIONAL INPUTS:
;	columns - vector of column numbers to be printed or a string
;		with column names separated by commas. If not supplied
;		or set to the null string, all columns are printed.
;
;	row1 - first row to print.  (default=0)
;	row2 - last row to print.  (default=last row in table)
;
; SIDE EFFECTS:
;	text is printed as directed by !textout
;
; HISTORY:
;	version 1, D. Lindler  Apr 89
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-----------------------------------------------------------------------------
	tab_read,name,tcb,tab		;read the table
	case n_params(0) of
		0: begin
			print,'TABLE_PRINT-- no table name suppled'
			retall
		   end
		1: tab_print,tcb,tab
		2: tab_print,tcb,tab,columns
		3: tab_print,tcb,tab,columns,row1
		4: tab_print,tcb,tab,columns,row1,row2
	endcase
	textclose
	return
	end
