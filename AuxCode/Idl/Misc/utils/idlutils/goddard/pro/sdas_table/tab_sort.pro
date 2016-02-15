pro tab_sort,column,tcb,tab
;+
; NAME:
;	TAB_SORT
; PURPOSE:
;	Procedure to sort table by the specified column
;
; CALLING SEQUENCE:
;	tab_sort, column, tcb, tab
;
; INPUTS:
;	column - column name or number to sort on
;	tcb - table control block
;
; INPUT/OUTPUTS:
;	tab - table array
;
; HISTORY:
;	version 1  D. Lindler  April 1989
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-----------------------------------------------------------------------
	values=tab_val(tcb,tab,column)	;extract column
	index=sort(values)		;sort column
	nrows=n_elements(index)		;number of rows
	newtab=tab
	for i=0,nrows-1 do newtab[0,i]=tab[*,index[i]]	;sort table
	tab=newtab
	return
	end
