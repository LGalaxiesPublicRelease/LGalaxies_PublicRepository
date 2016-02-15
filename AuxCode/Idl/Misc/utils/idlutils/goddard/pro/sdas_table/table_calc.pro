pro table_calc,table,expression,table_out
;+
; NAME:
;	TABLE_CALC 
; PURPOSE:
;	Adds a new table column from a expression using existing columns
;
; CALLING SEQUENCE:
;	table_calc, table, expression, table_out
;
; INPUTS:
;	table - input SDAS table table
;	expression - expression for new or updated column values.
;		Any legal IDL expression is valid where existing
;		column names can be used as variables.  User functions
;		within the expression are allowed if the function
;		is in an IDL library or previously compiled.
;
; OPTIONAL INPUT:
;	table_out - output table name.  If not supplied, the
;		input name is used.
;
; OUTPUTS:
;	a new SDAS table file is created.
;
; EXAMPLES:
;
;	 create a column WAVELENGTH in table TAB which is the average
;	of the WLOW and WHIGH columns of table TAB.
;
;		table_calc,'tab','WAVELENGTH=(WLOW+WHIGH)/2.0'
;
;	add a column SINX which is the sin of column X to table JUNK.
;
;		table_calc,'junk','SINX=sin(X)'
;
;	add 10.0 to an existing column in table MYTAB.
;
;		table_calc,'mytab','flux=flux+10.0'
;
; HISTORY
;	version 1  D. Lindler November, 1989
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;----------------------------------------------------------------------
;
; set output file name
;
	if n_params(0) lt 3 then table_out=table
;
; read input table
;
	tab_read,table,tcb,tab,h
;
; get name of the output column
;
	st = expression
	outcol = gettok(st,'=')
	if strtrim(st) eq '' then begin
		print,'TABLE_CALC--invalid expression supplied'
		retall
	endif
	st=strupcase(st)
;
; get names of the input columns and read if needed into a variable
; with a name equal to the column name.
;
	tab_size,tcb,nrows,ncols
	for i=0,ncols-1 do begin
	    tab_col,tcb,i+1,offset,width,dtype,name
	    name = strtrim(name)
	    if strpos(st,strupcase(name)) ge 0 then begin
		istat = execute(name+"=tab_val(tcb,tab,'"+name+"')")
		if istat le 0 then return
	    endif
	endfor
;
; execute the expression with the result going into varible x
;
	istat = execute('x='+st)
	if istat le 0 then begin
		print,'TABLE_CALC - error processing the expression'
		retall
	endif
;
; add column to the table
;
	tab_put,outcol,x,tcb,tab
;
; write the new table
;
	sxaddhist,' TABLE_CALC:  '+outcol+' = '+st,h
	tab_write,table_out,tcb,tab,h
	return
	end
