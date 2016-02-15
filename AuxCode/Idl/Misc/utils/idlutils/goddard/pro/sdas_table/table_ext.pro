pro table_ext,name,columns,v1,v2,v3,v4,v5,v6,v7,v8,v9
;+
; NAME:
;	TABLE_EXT
; PURPOSE:
;	Routine to extract columns from an STSDAS table
;
; CALLING SEQUENCE:
;	TABLE_EXT, name, columns, v1, [v2,v3,v4,v5,v6,v7,v8,v9]
; INPUTS:
;	name - table name, scalar string
;	columns - table columns to extract.  Can be either 
;		(1) String with names separated by commas
;		(2) Scalar or vector of column numbers
;
; OUTPUTS:
;	v1,...,v9 - values for the columns
;
; EXAMPLES:
;	Read wavelength and throughput vectors from STSDAS table, wfpc_f725.tab
;
;	IDL> table_ext,'wfpc_f725.tab','wavelength,throughput',w,t
;		or
;	IDL> table_ext,'wfpc_f725.tab',[1,2],w,t
;	
; PROCEDURES CALLED:
;	GETTOK(), TAB_READ, TAB_VAL()
; HISTORY:
;	version 1  D. Lindler  May 1989
;	Accept Column Numbers as well as names, W. Landsman  February 1996
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;---------------------------------------------------------------------
 if N_params() LT 3 then begin
	print,'Syntax -	TABLE_EXT, name, columns, v1, [v2,v3,v4,v5,v6,v7,v8,v9]'
	return
 endif
 n_ext = N_params() - 2
 sz_name = size(columns)
 strng = sz_name[sz_name[0]+1] EQ 7    ;Is colname a string?
 tab_read,name,tcb,tab
 colnames= columns
 for i = 0,n_ext-1 do begin
	if strng then colname = strtrim(gettok(colnames,','),2) $
		 else colname = colnames[i]
	v = tab_val(tcb,tab,colname)
	command = 'v'+strtrim(i+1,2)+'=v'
	istat = execute(command)
	endfor
 return
 end 


