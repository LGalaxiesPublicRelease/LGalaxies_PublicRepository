pro tab_write,name,tcb,tab,header
;+
; NAME:
;	TAB_WRITE
; PURPOSE:
;	Routine to write an stsdas table to disk
;
; CALLING SEQUENCE:
;	tab_write, name, tcb, tab, header
;
; INPUTS:
;	name - file name (default extension = .tab)
;	tcb - table control block
;	tab - table array
;
; OPTIONAL INPUT:
;	header - FITS header array
;
; HISTORY:
;	version 1  D. Lindler   April 1989
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;------------------------------------------------------------------------
;
; determine file size in blocks
;
	maxpar = (n_elements(header)-1)>0
	maxcol = tcb[6]
	ncols = tcb[5]
	rowlen = tcb[7]
	max_rowlen = tcb[8]>rowlen
	nrows = tcb[3]
	max_rows = tcb[4]>nrows
	nbytes = 12*4 + maxpar*80 + 16*4*maxcol + max_rows*max_rowlen*2
	nrecs = (nbytes+511)/512
;
; open output file
;
	file=strtrim(name,2)
	if strpos(file,'.') eq -1 then file=file+'.tab'
	if !version.os eq 'vms' then begin
		openw,unit,file,512,/none,/block,/get_lun
	endif else openw,unit,file,/get_lun
;
; write user parameters
;
	npar = 0
	if maxpar gt 0 then begin
		hblock = bytarr(80,maxpar)
		hblock[0:7,*] = 32b		;blanks
		for i = 0,maxpar-1 do begin
			line = header[i]
			keyword = strtrim(strmid(line,0,8))
			if keyword eq 'END' then goto,done_par
			if (keyword ne 'HISTORY') then $ 
				value = sxpar(header,keyword) $
				else value=strtrim(strmid(line,8,72))
			s = size(value) & valuetype = s[s[0]+1]
			case valuetype of
			  7: begin 			;string
				value = string(value)
				type = 't'
			     end
			  1: begin
				value = strtrim(value,2)	;byte / boolean
				type = 'b'
			     end
			  2: begin			;integer
				value = strtrim(value,2)
				type = 'i'
			     end
			  4: begin			;real*4
				value = strtrim(string(value,'(G16.8)'),2)
				type = 'r'
			     end
			  3: begin			;lonword integer
				value = strtrim(value,2)
				type = 'i'
			      end
			  5: begin			;real*8
				value = strtrim(string(value,'(G24.16)'),2)
				type = 'd'
			      end
			endcase
			hblock[0,i] = byte(keyword)
			hblock[8,i] = byte(type)
			hblock[9,i] = byte(value)
			npar = npar+1
		endfor	
done_par:	
		brec = assoc(unit,hblock,12*4)	
		brec[0] = hblock
	endif
;
; write column information
;
	if ncols gt 0 then begin
		brec = assoc(unit,tcb[*,1:ncols],12*4+maxpar*80)
		brec[0] = tcb[*,1:ncols]
	end
;
; write table
;
	if nrows*ncols gt 0 then begin
	    tab_pos = 12*4+maxpar*80+maxcol*64	;position of table within file
	    if tcb[9] ne 12 then begin		;row ordered table
		brec = assoc(unit,tab,tab_pos)
		brec[0] = tab
		tcb[9] = 11
	      end else begin			;column ordered table
		for i = 1,ncols do begin
			tab_col,tcb,i,offset,width
			position = offset*max_rows+tab_pos
			v = tab[offset:offset+width-1,*]
			brec = assoc(unit,v,position)
			brec[0] = v
		endfor
	    endelse
	endif
;
; write header
;
	brec = assoc(unit,bytarr(9*4))
	brec[0] = byte([npar,maxpar,nrows,max_rows,ncols, $
			maxcol,rowlen,max_rowlen,tcb[9]],0,9*4)
	close,unit
	free_lun,unit
return
end
