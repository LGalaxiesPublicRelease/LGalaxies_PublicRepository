pro table_help,tcb,header
;+
; NAME:
;	TABLE_HELP
; PURPOSE:
;	Procedure to decribe an SDAS table file.
;
; CALLING SEQUENCE:
;	table_help, tcb, header
;	table_help, name
;
; INPUTS:
;	tcb - table control block returned by TAB_READ or TAB_CREATE
;	name -	the table name
;
; OPTIONAL INPUTS:
;	header - header array returned by TAB_READ.  If supplied
;		it will be printed, otherwise it won't.
;
; SIDE EFFECTS:
;	text output as specified by !textout
;
; HISTORY:
;	version 1  D. Lindler  JAN 1988
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;----------------------------------------------------------------------
s=size(tcb) & type=s[s[0]+1]
if type eq 7 then begin		;file name supplied
	tab_read,tcb,tcbx,tab,header
   end else begin		;table control block supplied
	tcbx=tcb
end
;
; open output
;
textopen,'TABLE_HELP'
;
; print header info
;
ncol=tcbx[5,0]
printf,!textunit,'	'+string(byte(tcbx[*,tcbx[6,0]+1],0,64))
printf,!textunit,'         nrows=',tcbx[3,0],'               ncols=',tcbx[5,0]
printf,!textunit,'       maxrows=',tcbx[4,0],'             maxcols=',tcbx[6,0]
printf,!textunit,'          npar=',tcbx[1,0],'          row length=',tcbx[7,0],' words'
printf,!textunit,'     max n_par=',tcbx[2,0],'      max row length=',tcbx[8,0],' words'
ttype='Row-ordered'
if tcbx[9,0] eq 12 then ttype='Column-ordered'
printf,!textunit,'    Table type=  ',ttype
printf,!textunit,' '
;
; print column information
;
if ncol gt 0 then begin
   printf,!textunit,$
      ' col. number    column_name            units     type    SPP format code'
   printf,!textunit,' '
   for i=1,ncol do begin
	desc=tcbx[*,i]		;column description
	name=nulltrim(string(byte(desc,16,20)))
	units=nulltrim(string(byte(desc,36,20)))
	form=nulltrim(string(byte(desc,56,8)))
	if desc[3] lt 0 then dtype = 2 else dtype = desc[3]
	case dtype of
		6: type='Real*4'
		7: type='Real*8'
		4: type='Integer*4'
		1: type='Boolean'
		2: type='String'
	endcase
	printf,!textunit,'$(i6,a20,a20,a10,a10)',i,name,units,type,form
   end
end
if n_elements(header) gt 0 then begin
	printf,!textunit,' '
	for i=0,n_elements(header)-1 do begin
		printf,!textunit,header[i]
		if strmid(header[i],0,3) eq 'END' then goto,at_the_end
	endfor
at_the_end:
endif
textclose
return
end
