function tab_val, tcb, table, column, rows
;+
; NAME:
;	TAB_VAL
; PURPOSE:
;	Routine to read a column from an SDAS table file
;
; CALLING SEQUENCE:
;	values = tab_val( tcb, table, column, [ rows ] )
; INPUTS:
;	tcb - table control block returned by tab_val
;	table - table array returned by tab_val
;	column - scalar column name or number
; OPTIONAL INPUT:
;	rows - scalar giving row number or vector giving rows.
;		If not supplied all rows are returned.
; OUTPUT:
;	the values for the specified column (and rows) is returned
;	as the function value.  If row is specified as a scalar
;	(single row) then the result will be a scalar.
; HISTORY:
;	version 1  D. Lindler  Jan. 1988
;       Allow for a null column Landsman/Feggans    April 1992
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-------------------------------------------------------------------------
on_error,1         ;Return to caller
;
; get column information
;
tab_col,tcb,column,offset,width,datatype,name,units,format
if !err lt 0 then $
 	message,'Specified column, '+string(column) + ', not found'
;
nrows=tcb[3,0]
;
; determine row range
;
	if n_elements(rows) eq 0 then begin
		row1 = 0
		row2 = nrows-1
	   end else begin
		s = size(rows)
		if s[0] gt 0 then begin
			row1 = min(rows)
			row2 = max(rows)
		   end else begin
			row1 = rows
			row2 = rows
		end
	end
	nrows = row2-row1+1
;
; extract column and covert to correct data type
;
col=table[offset:offset+width-1,row1:row2]
case datatype of
	6: col=float(col,0,nrows)
	7: col=double(col,0,nrows)
	4: col=long(col,0,nrows)
	1: col=long(col,0,nrows)
	2: begin
		col=string(col)
		if nrows eq 1 then begin
			c = col
                        if strlen(col) gt 0 then $
				col=strarr(strlen(col),1) $
			else col = strarr(1,1) + string(32b)    ;Null field
			col[0]=c
		endif
		for i=0,nrows-1 do col[i]=nulltrim(col[i])
	   end
	else: message,'Unsupported column data type'
endcase
if n_elements(rows) gt 0 then return,col[rows-row1] else return,col
end
