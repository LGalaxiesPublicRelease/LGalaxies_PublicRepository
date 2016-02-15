pro sxhedit,name,h
;+
; NAME:
;	SXHEDIT                             
; PURPOSE:
;	Routine to interactively edit an STSDAS header on disk.
; EXPLANATION:
;	VMS: uses EDT.
;	Unix: uses whatever your EDITOR environment variable is set to.
;
; CALLING SEQUENCE:
;	sxhedit, name, [ h ]
;
; INPUTS:
;	name - header file name (default extension is .hhh)
;
; OUTPUTS:
;	h - (optional) edited header
;
; SIDE EFFECTS:
;	A new version of the file will be created.
;
; HISTORY:
;	Version 1  D. Lindler July  1987
;	Version 2  JAH Dec '88:  Converted to Sun IDL.
;	Modified   D. Neill Sept, 1990: Now deletes all versions of sxhedit.tmp
;			made compatable with Unix
;	Modified   D. Neill Apr, 1991: Ensures 80 char headers and will not
;		create new version if no changes made.
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;---------------------------------------------------------------------------
 if N_params() EQ 0 then begin
	print,'Syntax  - SXHEDIT, name, [ hdr ]
	return
 endif

; Check os

 do_vms=(!version.os eq 'vms')

; Read input header

 sxhread,name,h

; write to standard text file

 openw,unit,'sxhedit.tmp',/get_lun
 for i=0,n_elements(h)-1 do printf,unit,h[i]
 close,unit

; spawn to EDITOR

if do_vms then edt,'sxhedit.tmp' $
	  else if (getenv('EDITOR') eq '') then spawn,'vi sxhedit.tmp' $
					   else spawn,'$EDITOR sxhedit.tmp'

; read edited file

 he = strarr(n_elements(h)+100)
 openr,unit,'sxhedit.tmp',/delete
 st=''
 n=0
 while not eof(unit) do begin
	readf,unit,st

; ensure 80 characters in header lines

	if strlen(st) lt 80 then begin
		bl=string(replicate(32b,80))
		strput,bl,st
		st=bl
	endif else if strlen(st) gt 80 then st=strmid(st,0,80)
	he[n]=st
	n=n+1
 end
 he = he[0:n-1]
 free_lun,unit

; write new header (if needed)

 t = where((he ne h), dif)
 if dif gt 0 then begin
	sxhwrite,name,he
	sxhread,name,h

; delete original temp file

	if do_vms then spawn,'delete sxhedit.tmp;'
 endif

 return
 end
