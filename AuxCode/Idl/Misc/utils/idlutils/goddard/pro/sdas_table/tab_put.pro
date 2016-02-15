pro tab_put,column,values,tcb,tab,row
;+
; NAME:
;	TAB_PUT   
; PURPOSE:
;	Procedure to place new values into a STSDAS table.
;
; CALLING SEQUENCE:
;	tab_put, column, values, tcb, tab, row
;
; INPUTS:
;	column - column name or number (if it is a new column then
;		a column name must be specified)
;	values - data values to add to the table
;
; INPUT/OUTPUTS:
;	tcb - table control block
;	tab - table array
;
; OPTIONAL INPUT:
;	row - starting row to insert values
;
; HISTORY:
;	version 1  D. Lindler   April 1989
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;--------------------------------------------------------------------------
	if n_params(0) lt 5 then row = 0
;
; add column if it does not already exist
;
	tab_addcol,column,values,tcb,tab
;
; get column information
;
	tab_col,tcb,column,offset,width,datatype
	colnum = !err
;
; expand number of rows if the new data doesn't fit
;
	nrows = tcb[3]		;present number of rows
	nv = n_elements(values)
	last_row = nv+row-1
	if (last_row+1) gt nrows then begin
		tab_expand,tcb,tab,0,last_row+1
		tcb[3]=last_row+1		   ;new number of rows
		tab_nullrow,tcb,tab,nrows,last_row ;fill new rows with nulls
	endif
;
; convert data to correct type
;
	case datatype of
		1: v = long(values ne 0)
		2: v = string(values)
		4: v = long(values)
		6: v = float(values)
		7: v = double(values)
	endcase
	if datatype ne 2 then begin
		v = byte(v,0,width,nv)
	   end else begin
		v = byte(v)
		s = size(v) & ns=s[1]
		if ns gt width then v = v[0:width-1,*]
	endelse
;
; insert into table
;
	tab[offset,row] = v
return
end

