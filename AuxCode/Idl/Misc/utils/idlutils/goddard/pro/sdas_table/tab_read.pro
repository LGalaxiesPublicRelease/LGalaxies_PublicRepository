pro tab_read,name,tcb,table,header
;+
; NAME:
;   TAB_READ   
; PURPOSE:
;   Procedure to read an SDAS table file
; CALLING SEQUENCE:
;	tab_read,name,tcb,table,header
; INPUTS:
;	name - name of the table file
; OUTPUTS:
;	tcb - table control block 
;		Longword array of size 16 x maxcols+2
;		where maxcols is the maximum number of columns
;		allocated for the table.
;		tcb(*,0) contains:
;		   word	0	SPARE
;			1	number of user parameters
;			2	max. number of user par. allowed
;			3	number of rows in the table
;			4	number of allocated rows (for col. ordered tab)
;			5	number of columns defined
;			6	max number of columns
;			7	length of row used (in units of 2-bytes)
;			8	max row length (in units of 2-bytes)
;					relevant only for row ordered tables.
;			9	table type (11 for row order, 12 for col. order)
;			15	update flag (0-readonly, 1-update)
;		tcb(*,i) contains description of column i
;		   word 0	column number
;			1	offset for start of row in units of 2-bytes
;			2	width or column in 2-byte units
;			3	data type
;					6 = real*4
;					7 = real*8
;					4 = integer*4
;					1 = boolean*4
;					2 = character string
;			4-8	ascii column name up to 19 characters
;			9-13	column units (up to 19 characters)
;			14-15	format string
;		tcb(*,max number of columns+1)= file name
;
;	table - table array, Byte array row length (bytes) x nrows
;
; 	header - header parameters in form usable by sxpar, sxaddhist,
;		sxaddpar, ect.
; HISTORY:
;	Version 1  D. Lindler  Jan 88
;	Converted to NEW IDL  April 90  D. Lindler
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;------------------------------------------------------------------------
if n_params(0) lt 2 then begin
	print,'TAB_READ -- output table control block not supplied'
	retall
endif
;
; open file
;
fdecomp,name,disk,dir,root,ext
if ext eq '' then ext = 'tab'
fname=disk+dir+root+'.'+ext
!err=0
openr,unit,fname,/block,/get_lun
if !err lt 0 then begin
	free_lun,unit
	message,'Error opening file '+fname,/CON
	print,!ERROR_STATE.MSG
	retall
end
;
; read header record
;
lrec=assoc(unit,lonarr(12),0)
h=lrec[0]
maxpar=h[1]
maxcol=h[5]
type=h[8]
;
; create table control block
;
tcb=lonarr(16,maxcol+2)
tcb[0]=[unit,h,0,0,0]
sbyte=12*4+h[1]*80L			;starting byte of column descriptions
lrec = assoc(unit,lonarr(16,maxcol),sbyte)
tcb[0,1]=lrec[0]			;read col. descriptions
byte_name=replicate(32b,64)
byte_name[0]=byte(fname)
tcb[0,maxcol+1]=long(byte_name,0,16)	;place filename into tcb
;
; read data
;
rowlen=MAX([tcb[7,0],tcb[8,0]])*2
nrows=MAX([tcb[4,0],tcb[3,0]])
offset=12*4+tcb[2,0]*80+tcb[6,0]*16*4	;offset of data in the file

brec=assoc(unit,bytarr(1))		; we will read into a byte array

if type eq 11 then begin	; row ordered
	brec = assoc(unit,bytarr(rowlen,nrows,/nozero),offset)
	table=brec[0]		;read table
   end else begin		; column ordered table
	table=bytarr(rowlen,nrows)
	ncols=tcb[5,0]
	if ncols gt 0 then begin
		for i=1,ncols do begin
			tab_col,tcb,i,position,width,datatype
			brec = assoc(unit,bytarr(width,nrows,/nozero),offset)
			column=brec[0]			 ; read column
			offset=offset+nrows*width	;offset to next col
			table[position,0]=column	;insert into table
		endfor
	endif
endelse
;
; read header
;
if n_params(0) gt 3 then begin	;read only if output parameter supplied
	npar=tcb[1,0]
	header=strarr(npar+1)
	header[0]='END'+string(replicate(32b,77))
	if npar gt 0 then begin
	    brec = assoc(unit,bytarr(80,npar,/nozero),12*4)
	    htab=string(brec[0])		;read table header
	    for i=0,npar-1 do begin		;convert to FITS header
		line=htab[i]
		keyword=strmid(line,0,8)
		type=strmid(line,8,1)
		value=nulltrim(strmid(line,9,71))
		if strtrim(keyword ne '') then begin
		  case type of			;convert value to correct type
			't':
			'b': if value eq '0' then value='F' else value='T'
			'i': value=long(strtrim(value))
			'r': value=float(strtrim(value))
			'd': value=double(strtrim(value))
		  endcase
		  sxaddpar,header,keyword,value
		end
	    endfor
	endif
endif
close,unit
free_lun,unit
return
end
