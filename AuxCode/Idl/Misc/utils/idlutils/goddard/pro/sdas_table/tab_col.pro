pro tab_col,tcb,column,offset,width,datatype,name,units,format
;+
; NAME:
;	TAB_COL  
; PURPOSE:
;	Procedure to extract column information from table control block
;
; CALLING SEQUENCE:
;	tab_col, tcb, column, offset, width, datatype, name, units, format
;
; INPUTS:
;	tcb - table control block returned by tab_open.
;	column - column name (string) or column number
;
; OUTPUTS:
;	offset - column offset bytes
;	width - column width in bytes
;	datatype - column data type:
;		6 - real*4
;		7 - real*8
;		4 - integer*4
;		1 - boolean
;		2 - character string
;	name - column name
;	units - column units
;	format - format code
;
; SIDE EFFECTS:
;	If the column is not found then !err is set to -1.
;	Otherwise !err is set to the column number (starting at one).
;
; HISTORY:
;	version 1  D. Lindler  Jan 88
;	Converted to NEW IDL  April 90
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;----------------------------------------------------------------------------
;
; determine if valid control block
;
s=size(tcb)
ndim=s[0]
if (ndim ne 2) or (s[1] ne 16) or (s[ndim+1] ne 3) then begin
	Print,'TAB_COL -- invalid table control block'
	print,'It must be a 2-D long word array with first dimension=16'
	retall
endif
;
; get number of columns in the table
;
ncols=tcb[5,0]
;
; determine if column name of number supplied
;
s=size(column)
ndim=s[0]
if ndim ne 0 then begin
	print,'TAB_COL -- column must be a scalar string or number'
	retall
endif
if s[ndim+1] ne 7 then begin	;number supplied
	colnum=long(column)
	if (colnum lt 1) or (colnum gt ncols) then begin
		print,'TAB_COL -- Invalid column number specified'
		print,'It must be between   1 and',ncols
		!err=-1
		return
	endif
    end else begin		;name specified
;
; loop and find name in control block
;
	cname=strupcase(strtrim(column))
	for i=1,ncols do begin
		name=nulltrim(string(byte(tcb[4:8,i],0,19)))
		if(strupcase(name) eq cname)then goto,found
	endfor
	!err=-1			;not found
	return
found:
	colnum=i
end
;
; extract information
;
offset=tcb[1,colnum]*2
width=tcb[2,colnum]*2
datatype=tcb[3,colnum]
if datatype lt 0 then begin
	width = -datatype
	datatype = 2
end
name=nulltrim(string(byte(tcb[4:8,colnum],0,19)))
units=nulltrim(string(byte(tcb[9:13,colnum],0,19)))
format=nulltrim(string(byte(tcb[14:15,colnum],0,8)))
!err=colnum
return
end	
