pro tab_print,tcb,tab,columns,row1,row2
;+
; NAME:
;	TAB_PRINT 
; PURPOSE:
;	Routine to print an stsdas table.
;
; CALLING SEQUENCE:
;	tab_print, tcb, tab, columns, row1, row2
;
; INPUTS:
;	tcb - table control block returned by TAB_READ
;	tab - table array read by TAB_READ
;
; OPTIONAL INPUTS:
;	columns - vector of column numbers to be printed or a string
;		with column names separated by commas. If not supplied
;		or set to the null string, all columns are printed.
;
;	row1 - first row to print.  (default=0)
;	row2 - last row to print.  (default=last row in table)
;
; SIDE EFFECTS:
;	text is printed as directed by !textout
;
; HISTORY:
;	version 1, D. Lindler  Apr 89
;	April 90  Converted to NEW IDL D. Lindler
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-----------------------------------------------------------------------------
;
; get table size
;
	nrows=tcb[3]
	ncols=tcb[5]
	maxcols=tcb[6]
	if nrows*ncols eq 0 then begin
		print,'TAB_PRINT--input table is empty'
		return
	endif
;
; set defaults for non-supplied parameters
;
	if n_params(0) lt 3 then columns=''
	if n_params(0) lt 4 then row1=0
	if n_params(0) lt 5 then row2=nrows-1
;
; convert columns to integer vector if not supplied that way
;
	s=size(columns) & dtype=s[s[0]+1]
	if dtype eq 7 then begin
		st=strtrim(columns)
		if st eq '' then begin
			colno=indgen(ncols)+1
		  end else begin
			colno=intarr(100)
			ncol=0
			while st ne '' do begin
			    next_col=strtrim(gettok(st,','),2)
			    tab_col,tcb,next_col
			    if !err lt 1 then begin
				print,'TAB_PRINT-- column '+next_col+ $
						' not found'
				retall
			    endif
			    colno[ncol]=!err
			    ncol=ncol+1
			endwhile
			colno=colno[0:ncol-1]
		endelse
	   end else begin
		colno=intarr(n_elements(columns))+columns
	endelse
	if (max(colno) gt ncols) or (min(colno) lt 1) then begin
		print,'TAB_PRINT-- Invalid column number specified'
		retall
	endif
;
; check validity of row numbers
;
	row1=row1>0
	row2=row2>row1<(nrows-1)
;
; get column information
;
	ncol=n_elements(colno)
	byte1=intarr(ncol)		;byte positions within row
	byte2=intarr(ncol)
	dtypes=intarr(ncol)
	names=strarr(ncol)
	units=strarr(ncol)
	formats=strarr(ncol)		;fortran formats
	fwidth=intarr(ncol)		;print field widths
	for i=0,ncol-1 do begin
		tab_col,tcb,colno[i],off,width,dtype,name,unit,sppformat
		byte1[i]=off
		byte2[i]=off+width-1
		dtypes[i]=dtype
		names[i]=name
		units[i]=unit
		tab_spptofor,sppformat,format,width	;convert to fortran
		fwidth[i]=width
		formats[i]='('+format+')'
	endfor
;
; create header lines
;
	line1='  row  '
	line2='       '
	for i=0,ncol-1 do begin
		name=strmid(strtrim(names[i]),0,fwidth[i])
		while strlen(name) lt fwidth[i] do name=' '+name
		line1=line1+' '+name
		unit=strmid(strtrim(units[i]),0,fwidth[i])
		while strlen(unit) lt fwidth[i] do unit=' '+unit
		line2=line2+' '+unit
	endfor
	line0='    '+nulltrim(string(byte(tcb[*,maxcols+1],0,64)));file name
;
; loop on rows and print headers every 50 lines
;
	textopen,'tab_print'
	nlines=50
	for row=row1,row2 do begin

		if nlines ge 50 then begin	;new page ?
			printf,!textunit,string(12B)	;page eject
			printf,!textunit,line0
			printf,!textunit,line1
			printf,!textunit,line2		
			printf,!textunit,' '
			nlines=0
		endif
		
		line=string(row,'(i7)')
		for i=0,ncol-1 do begin
			val=tab[byte1[i]:byte2[i],row]
			case dtypes[i] of
				6: val=float(val,0)
				7: val=double(val,0)
				4: val=long(val,0)
				1: val=long(val,0)
				2: val=nulltrim(string(val))
			endcase
			line=line+' '+string(val,formats[i])
		endfor
		printf,!textunit,line
		nlines=nlines+1
	endfor
	textclose
	return
	end
