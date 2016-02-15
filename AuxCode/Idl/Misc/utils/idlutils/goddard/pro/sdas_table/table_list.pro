pro table_list,name,row1,row2,header=header,textout=textout
;+
; NAME:
;	TABLE_LIST  
; PURPOSE:
;	List the contents of an STSDAS table.
; EXPLANATION:
;	Procedure to list contents of an STSDAS table.  This does not
;	print the table in tabular form but instead for each row
;	prints the column name followed by its value (one column per
;	output line.
;
; CALLING SEQUENCE:
;	table_list, name, row1, row2, [ TEXTOUT=, /HEADER ]
;
; INPUTS:
;	name - table name
;
; OPTIONAL KEYWORD INPUT:
;	TEXTOUT  - Scalar string giving output file name, or integer (1-5)
;		specifying output device.   See TEXTOPEN for more info.
;		Default is to display output at the terminal
;	HEADER - if set, the header is printed before the selected row printout
;
; OPTIONAL INPUTS:
;	row1 - first row to list (default = first row)
;	row2 - last row to list (default = last row)
;
; OUTPUT:
;	text output is written to the output device specified by the TEXTOUT
;	keyword, or the nonstandard system variable !TEXTOUT
;
; PROCEDURES USED:
;	TAB_COL, TAB_READ, TAB_SIZE, TAB_VAL(), TEXTOPEN, TEXTCLOSE
;
; HISTORY:
;	version 1  D. Lindler   May 1989
;	July 1996, DJL, added /header keyword to optionally print header
;	August 1996, WBL, added TEXTOUT keyword
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;----------------------------------------------------------------------------
	if N_params() LT 1 then begin
	print,'Syntax = table_list, name, [ row1, row2, TEXTOUT= ,/HEADER]'
	return
	endif
     
	tab_read,name,tcb,tab,h			;read table
	tab_size,tcb,nrows,ncols		;get its size

	if n_params(0) lt 2 then row1=0		;determine rows to list
	if n_params(0) lt 3 then row2=nrows-1
	if not keyword_set(TEXTOUT) then textout = !TEXTOUT

	row2=row2>row1

	textopen,'table_list',textout=textout		;open text file
	if keyword_set(header) then printf,!textunit,h

	for i=row1,row2 do begin
	    printf,!textunit,'--------------------------------------- Row ',i
	    for j=1,ncols do begin
		tab_col,tcb,j,off,width,dtype,colname	;get column name
		value=tab_val(tcb,tab,j,i)		;get its value
		printf,!textunit,string(colname,'(A25)')+' = '+string(value)
	    end
	    printf,!textunit,' '
	end
	textclose,textout=textout
return
end
