pro tab_addcol,name,data,tcb,tab
;+
; NAME:
;	TAB_ADDCOL   
; PURPOSE:
;	Procedure to add a new column to an existing STSDAS table.
;
; CALLING SEQUENCE:
;	tab_addcol, name, data, tcb, tab
;
; INPUTS:
;	name - column name
;	data - sample data of type to be written to the column.
;		This parameter is only used to determine data type.
;
; INPUT/OUTPUTS:
;	tcb - table control block
;	tab - table array
;
; HISTORY:
;	version 1  D. Lindler April 89
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;------------------------------------------------------------------------
;
; If column alread exists, do not add it
;
	tab_col,tcb,name
	if !err gt 0 then return
;
; determine data type and change to stsdas table code
;
	s=size(data) & dtype=s[s[0]+1] & ndim=s[0]
	case dtype of
		0: type=6		;undefined, assume real*4
		7: type=2		;string
		1: type=1		;byte (assume boolean)
		2: type=4		;integer
		4: type=6		;real*4
		3: type=4		;integer
		5: type=7		;double precision
		6: begin
			print,'TAB_ADDCOL-- complex data type not supported'
			retall
		    end
	endcase
;
; determine field width
;
	case type of
		2: begin		;string
			width=max(strlen(data))
			width=(width+3)/4*2 ;change to multiple of 4 bytes
		   end
		7: width=4		;double precision
		else: width=2		;integer, real, or boolean
	endcase
;
; get current table size
;
	nrows=tcb[3]
	ncols=tcb[5]
	rowlen=tcb[7]
	maxrows=tcb[4]>nrows
	maxcols=tcb[5]
	maxrowlen=tcb[8]>rowlen
;
; updated table size
;
	new_nrows=nrows>n_elements(data)
	new_ncols=ncols+1
	new_rowlen=rowlen+width
;
; do we need to expand the table
;
	if (new_nrows gt maxrows) or (new_ncols gt maxcols) or $
	   (new_rowlen gt maxrowlen) then $
		tab_expand,tcb,tab,new_ncols,new_nrows,new_rowlen
;
; construct default print format
;
	case type of
		1: form='8b'				;boolean
		2: form=''+strtrim(width*2,2)+'s'	;string
		4: form='11d'				;integer
		6: form='16.8g'			;real*4
		7: form='20.12g'			;real*8
	endcase
;
; insert column information into tcb
;
	newcol=bytarr(48)
	newcol[0]=byte(name)
	newcol[40]=byte(form)
	tcb[0,new_ncols]=[long(new_ncols),rowlen,width,type,long(newcol,0,12)]
	tcb[5]=new_ncols
	tcb[7]=new_rowlen
;
; fill new column with nulls
;
	if nrows gt 0 then begin
		case type of
			1: null=0L
			2: null=string(0b)
			4: null='80000000'XL
			6: null=1.6e38
			7: null=1.6d38
		endcase
		if type eq 2 then nulls=bytarr(width*2,nrows)+null $
			     else nulls=replicate(null,nrows)
		nulls=byte(nulls,0,width*2,nrows)
		tab[rowlen*2,0]=nulls
	endif
return
end
