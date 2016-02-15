pro table_append,list,name
;+
; NAME:
;	TABLE_APPEND
; PURPOSE:
;	Routine to append STSDAS tables to create a single table. 
;	Input tables must all have identical columns.
;
; CALLING SEQUENCE:
;	table_append,list,name
;
; INPUTS:
;	list - string array listing the file names or a string
;		scalar giving a file name template.
;	name - output file name.
; SIDE EFFECTS:
;	a new STSDAS table is created with the specified name.
;
; OPERATIONAL NOTES:
;	all input tables must have the same number of columns
;	with the same names, datatypes, and column order.
;	Header parameters are taken only from the first table.
;
; HISTORY:
;	version 1  D. Lindler	April 1989
;       Removed call to non-standard system variable !DUMP WBL  September 1997
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;---------------------------------------------------------------------
;
; get list of file names if list is a scalar
;
	if n_elements(list) eq 1 then files=findfile(list) $
				 else files=list
		print,'Merging files:'
		print,files
	nfiles=n_elements(files)
;
; read first table as template for the output table
;
	tab_read,strtrim(files[0]),tcb,tab,h
	nrows=tcb[3]
	max_rows=tcb[4]>nrows
	ncols=tcb[5]
	rowlen=tcb[7]
	max_rowlen=tcb[8]>rowlen
	cnames=strarr(20,ncols)
	for i=1,ncols do begin
		tab_col,tcb,i,offset,width,datatype,cname
		cnames[i-1]=cname
	endfor
;
; loop on remaining tables
;
	if nfiles gt 1 then begin
	    for ifile=1,nfiles-1 do begin
		tab_read,strtrim(files[ifile]),tcb1,tab1
		nr=tcb1[3]			;number of rows
		nc=tcb1[5]			;number of columns
		rl=tcb1[7]			;row length
		if (rl ne rowlen) or (nc ne ncols) then begin
			print,'table_append-- input tables are not compatible'
			retall
		endif
;
; do columns match
;
		for i=1,ncols-1 do begin
			tab_col,tcb1,i,offset,width,datatype,cname
			if strtrim(cname) ne strtrim(cnames[i-1]) then begin
				print,'table_append-- column names must be'+ $
					' the same in all input tables'
				retall
			endif
		endfor
;
; do we need to expand the output table?
;
		if (nr+nrows) gt max_rows then $
			tab_expand,tcb,tab,0,nrows+nr*(nfiles-ifile+1)
		max_rows=tcb[4]
;
; insert new rows and update the table control block
;
		tab[0,nrows]=tab1[0:rowlen*2-1,0:nr-1]
		nrows=nrows+nr
		tcb[3]=nrows
	    endfor
	endif
;
; write the new table
;
	tab_write,name,tcb,tab,h
	return
end
