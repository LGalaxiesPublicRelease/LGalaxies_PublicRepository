pro tab_nullrow,tcb,tab,row1,row2
;+
; NAME:
;	TAB_NULLROW 
; PURPOSE:
;	Insert null row(s) into a STSDAS table
;
; CALLING SEQUENCE:
;	tab_nullrow, tcb, tab, [ row1, row2  ]
;
; INPUTS:
;	tcb - table control block
;
; INPUT/OUTPUTS:
;	tab - table array
;
; OPTIONAL INPUTS:
;	row1 - first row number to insert nulls (default=0)
;	row2 - last row number to insert nulls (default = last row)
;
; HISTORY:
;	version 1, D. Lindler  Apr 89
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;------------------------------------------------------------------------------
;
; get table size
;
	nrows=tcb[3]
	ncols=tcb[5]
	rowlen=tcb[7]
	if ncols eq 0 then return	;wasn't that easy.
;
; set defaults
;
	if n_params(0) lt 3 then row1=0
	if n_params(0) lt 4 then row2=nrows-1
	row1=row1>0
	row2=row2<(nrows-1)
	if row2 lt row1 then return	;that was easy too.
;
; create null row
;
	nullrow=bytarr(rowlen*2)
	for i=1,ncols do begin
		tab_col,tcb,i,offset,width,datatype
		if datatype gt 3 then begin
			case datatype of
				4: null='80000000'XL
				6: null=1.6E38
				7: null=1.6D38
			endcase
			nullrow[offset]=byte(null,0,width)
		endif
	endfor
;
; insert null rows
;
	for i=row1,row2 do tab[0,i]=nullrow
return
end
