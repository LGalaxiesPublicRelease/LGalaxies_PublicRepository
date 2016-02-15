pro sxhread, name, header
;+
; NAME:
;	SXHREAD                         
; PURPOSE:
;	Procedure to read a STSDAS header from disk.  
; EXPLANATION:
;	This version of SXHREAD can read three types of disk files
;	(1)  VMS Fixed record length 80 byte files, or GEIS files with
;		VMS buckets
;	(2)  Unix stream files with a CR after every 80 bytes
;	(3)  Variable length record files (Unix or VMS)
;
; CALLING SEQUENCE:
;	sxhread, name, header
;
; INPUT:
;	name - file name, scalar string.  An extension of .hhh is appended 
;		if not already supplied.   (Note STSDAS headers are required
;		to have a 3 letter extension ending in 'h'.)
; OUTPUT:
;	header - STSDAS header, string array
; NOTES:
;	SXHREAD  does not do any checking to see if the file is a valid
;	STSDAS header.    It simply reads the file into a string array with
;	80 byte elements
;
; HISTORY:
;	Version 1  D. Lindler  July, 1987
;	Version 2  M. Greason, August 1990
;	Use READU for certain ST VAX GEIS files   W. Landsman January, 1992
;	Read variable length Unix files  E. Deutsch/W. Landsman November, 1994
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;--------------------------------------------------------------------
 On_error,2                              ;Return to caller

 if N_params() LT 2 then $            
      message,'Syntax - SXHREAD, name, header',/NONAME

; Add extension name if needed

 hname = strtrim(name,2)
 if strpos(hname,'.',strpos(hname,']') ) EQ -1 then hname = hname + '.hhh'

 openr, unit, hname, /GET_LUN, ERROR = err    
 if err LT 0 then goto, BADFILE 

 len = 80  & ai = 99                    ;Usual header length is 80 bytes
 if !Version.os NE "vms" then begin     ;but Unix files may have an 
                                        ;embedded carriage returns to make
   atmp = assoc(unit,bytarr(85))           ;header length 81 bytes
   a=atmp[0] & ai=0
   while (a[ai] ne 10) and (a[ai] ne 13) and (ai lt 84) do ai=ai+1
   if (ai EQ 80) then len=81
   Point_lun, unit, 0            ;Back to the beginning of the file

 endif

; Get the number of lines in the header

 status = fstat(unit)
 nlines = status.size/len                      ;Number of lines in file
 if (!VERSION.OS EQ "vms") and (status.rec_len NE 80) then goto, VAR_LENGTH
 if (ai lt 80) then goto,VAR_LENGTH

; Read header

 header =  bytarr(len,nlines ,/NOZERO)
 On_ioerror, VAR_LENGTH        ;READU cannot be used on variable length records
 readu, unit, header
 header = string(header)
 On_ioerror,NULL 

 free_lun,unit             ;Close and free file unit

; Trim to the END line, and delete carriage returns if necessary

 endline = where( strmid(header,0,8) EQ 'END     ',nfound)
 if nfound gt 0 then header = header[0:endline[0]] else $
     message,'WARNING: No END statement found in header',/inform
 if len EQ 81 then header = strmid(header,0,80)
 return

VAR_LENGTH:                 ;Now try to read as variable length records

 Point_lun, unit, 0          ;Back to the beginning of file
 h = ''  & header = strarr( nlines)
 i = 0

 On_ioerror,NOEND            ;Can't use EOF function on certain GEIS files
 while ( strtrim( strmid(h,0,8), 2) NE 'END') do begin
    readf, unit, h
    if (strlen(h) LT 80) then h=h+string(replicate(32b,80-strlen(h)))
    header[i] = h                  ;Swapped with line above 95-Aug
    i = i + 1
    if i EQ nlines then begin 
            header = [header,strarr(100)]
            nlines = nlines + 100
     endif
 endwhile
 header = header[0:i-1]
 free_lun,unit
 return

NOEND: 
   message,'WARNING - No END statement found in header', /INFORM 
   free_lun,unit
   return

BADFILE: 
   if !VERSION.OS EQ 'vms' then hname = strupcase(hname)
   message,'Error opening file ' + ' ' + hname
   return

end
