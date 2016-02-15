Function OSFCNVRT,lname
;+
; NAME:
;	OSFCNVRT
;
; PURPOSE:
;	Return the correctly formatted logical directory syntax for the host OS
;
; CALLING SEQUENCE:
;	OSFCNVRT,lname
;
; INPUTS:
;	lname	- the file specification as a logical name + file name string
;
; OUTPUTS:
;	Returns appropriate string.
;
; SIDE EFFECTS:
;	None.
;
; RESTRICTIONS:
;	Assumes that the input is composed of only a logical and a filename combination
;	without lower directory garbage.
;
; PROCEDURE:
;	The operating system in !version.os is checked. If it equals:
;
;		'vms'		then a ':' is appended.
;
;		else		unix os is assumed and the logical portion is
;				uppercased, a '$' is prepended and a '/' is
;				appended.
;
; MODIFICATION HISTORY:
;	Written, JDNeill, May, 1990.
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;------------------------------------------------------------------------------
	name=strtrim(lname,2)
	if !version.os eq 'vms' then begin
		if (strpos(name,':') ge 0) or (strpos(name,'$') le 0) $
			then return,name
		if strpos(name,'$') eq 0 then $
			name=strmid(name,1,strlen(name)-2)
		log=gettok(name,'/')
		return,log+':'+name
	endif else begin
		if (strpos(name,'$') eq 0) or (strpos(name,':') le 0) $
			then return,name
		log=gettok(name,':')
		return,'$'+strupcase(log)+'/'+name
	endelse
;
	end	; osfcnvrt.pro
