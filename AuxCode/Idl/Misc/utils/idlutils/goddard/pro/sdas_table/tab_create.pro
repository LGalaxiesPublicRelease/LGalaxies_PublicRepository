pro tab_create,tcb,tab,maxcol,maxrows,row_len,tb_type
;+
; NAME:
;	TAB_CREATE  
; PURPOSE:
;	Procedure to create a new table file.
;
; CALLING SEQUENCE:
;	tab_create, tcb, tab, maxcol, maxrows, row_len, tb_type
;
; OUTPUTS:
;	tcb - table control block for reading from and writing
;		to the file (see tab_open for description)
;	tab - table array
;
; OPTIONAL INPUTS:
;	maxcol - maximum allocated number of columns [default=10]
;	maxrows - maximum allocated number of rows   [default=100]
;	row_len - row length in 2 byte units	     [default=2*maxcol]
;	tb_type - table type 'row' or 'column' ordered
;
; SIDE EFFECTS:
;	Table file is created and left opened to unit number tcb(0,0)
;	for writing.
;
; HISTORY:
;	version 1   D. Lindler   Dec. 88
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-------------------------------------------------------------------------
;
; set default parameters
;
	npar=n_params(0)
	if npar lt 3 then maxcol=10
	if npar lt 4 then maxrows=100
	if npar lt 5 then row_len=maxcol*2
	if npar lt 6 then tb_type='column'
	if strupcase(strmid(tb_type,0,1)) eq 'R' then typecode=11 else typecode=12
;
; create table control block
;
	tcb=lonarr(16,maxcol+2)
	tcb[0,0]=[0,0,0,0,maxrows,0,maxcol,0,row_len,typecode,0,0,0,0,0,1]
	byte_name=replicate(32b,64)
	tcb[0,maxcol+1]=long(byte_name,0,16)
	tab=bytarr(row_len*2,maxrows)
return
end
