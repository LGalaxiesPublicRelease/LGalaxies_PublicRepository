pro tab_expand,tcb,tab,maxcol,maxrow,rowlen
;+
; NAME:
;	TAB_EXPAND
; PURPOSE:
;	routine to expand the size of an SDAS table file.
;
; CALLING SEQUENCE:
;	tab_expand, tcb, tab, maxcol, maxrow, rowlen
;
; INPUT/OUTPUT:
;	tcb - table control block returned by routine TAB_READ
;		or TAB_CREATE.
;	tab - table array
;
; OPTIONAL INPUTS:
;	maxcol - new maximum number of columns.
;	maxrow - new maximum number of rows.
;	rowlen - new maximum row length in 2 byte units.
;
;	If maxcol, maxrow, or rowlen are supplied with
;	values less than the previous maximums, the previous
;	maximums are used.  All values are defaulted to zero
;	if not supplied.
;
; HISTORY:
;	Version 1   D. Lindler   Dec. 88
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-----------------------------------------------------------------------
;
; set default parameters
;
	npar=n_params(0)
	if npar lt 3 then return		;nothing to expand
	if npar lt 4 then maxrow=0
	if npar lt 5 then rowlen=0
;
; get old table sizes and parameters
;
	old_maxcol=tcb[6,0]
	old_maxrow=tcb[4,0]>tcb[3,0]
	old_rowlen=tcb[8,0]>tcb[7,0]
;
; set new table sizes
;
	new_maxcol=old_maxcol>maxcol
	new_maxrow=old_maxrow>maxrow
	new_rowlen=old_rowlen>rowlen
;
; increase size of table control block in maxcol increased
;
	if new_maxcol gt old_maxcol then begin
		new_tcb=lonarr(16,new_maxcol+2)
		new_tcb[0,0]=tcb
		new_tcb[0,new_maxcol+1]=tcb[*,old_maxcol+1]
		tcb=new_tcb
	endif
;
; increase table size if needed
;
	if (new_maxrow gt old_maxrow) or (new_rowlen gt old_rowlen) then begin
		oldtab=tab
		tab=bytarr(new_rowlen*2,new_maxrow)
		tab[0,0]=oldtab
	endif
;
; update control block
;
	tcb[6,0]=new_maxcol
	tcb[8,0]=new_rowlen
	tcb[4,0]=new_maxrow
	return
	end
