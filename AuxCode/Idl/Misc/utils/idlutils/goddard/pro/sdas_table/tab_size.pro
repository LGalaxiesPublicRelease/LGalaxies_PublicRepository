pro tab_size,tcb,nrows,ncols,maxrows,maxcols,rowlen,max_rowlen
;+
; NAME:
;	TAB_SIZE   
; PURPOSE:
;	Routine to extract the table size from a table control block
;
; CALLING SEQUENCE:
;	tab_size, tcb, nrows, ncols, maxrows, maxcols, rowlen, max_rowlen
;
; INPUTS:
;	tcb - table control block
;
; OUTPUTS:
;	nrows - number of rows in the table
;	ncols - number of columns in the table
;	maxrows - number of rows allocated
;	maxcols - number of columns allocated
;	rowlen - length of the rows in bytes
;	max_rowlen - allocated row length
;
; HISTORY:
;	version 1  D. Lindler  April 1989
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-------------------------------------------------------------------------
	nrows=tcb[3]
	maxrows=tcb[4]>nrows
	ncols=tcb[5]
	maxcols=tcb[6]
	rowlen=tcb[7]*2
	max_rowlen=(tcb[8]*2)>rowlen
	return
	end
