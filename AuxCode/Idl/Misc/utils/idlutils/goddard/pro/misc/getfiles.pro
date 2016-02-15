pro getfiles,list
;+
; NAME:
;	GETFILES
; PURPOSE:
;	Prompt the user to interactively specify a list of files
; EXPLANATION:
;	User can specify a single file per line or a range of files 
;	separated by a dash or comma.    Used, for example, by FITSRD to
;	return a list of file numbers on tape to read
;
; CALLING SEQUENCE:
;	getfiles, list
;
; OUTPUT:
;	LIST - integer array containing file numbers
;
; SIDE EFFFECTS:
;	User will be prompted to enter a list of file numbers
;
; REVISION HISTORY
;	Written D. Lindler November, 1985
;	Converted to Version 2 IDL,  August 1990
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
START:  
 list = intarr(2000)
 nlist = 0
 print,'Enter list of file numbers, one file or file range per line.'
 print,'Enter 0 to stop or H for help.  Type X to start over.'
 for i = 1,2000 do begin
  st = ''
  read,st
;
; check if aborted
;
  if strupcase(st) eq 'X' then goto,start
;
; check if help was asked for
;
 if (st eq 'h') or (st eq 'H') then begin
   print,' '
   print,'_____________________________________________________'
   print,' '
   print,'Files are entered one file number per line or'
   print,'   a range of files specified by the first and last'
   print,'   file number separated by a dash or comma.'
   print,'Type 0 (zero) to stop'
   print,'Type X to start over'
   print,' '
   print,'For example; if you type:
   print,'   5'
   print,'   8'
   print,'   12-18'
   print,'   0'
   print,' '
   print,' files 5, 8, and 12 through 18 would be specified'
   print,'  '
   print,'_____________________________________________________'
   print,' '
 end else begin
;
; check if a range was specified
;
   pos = strpos(st,'-')
   if pos lt 0 then pos = strpos(st,',')
   if pos ge 0 then begin		;range was specified
     len = strlen(st)			;get length of string
     first = fix(strmid(st,0,pos))	;get first file number
     last = fix(strmid(st,pos+1,len-pos-1));get last file number
     for file = first,last do begin
       list[nlist] = file
       nlist = nlist+1
     end; for file
;
; only single file specified
;
    end else begin
      file = string(st)
      if file eq 0 then goto,finished	;zero specified (done)
      list[nlist] = file
      nlist = nlist+1
    end; if pos ge 0
  end; if st eq 'H'
 end; for i

FINISHED:  
 list = list[0:nlist-1]

 return
 end
