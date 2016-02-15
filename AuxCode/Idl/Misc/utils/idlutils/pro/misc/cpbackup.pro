;+
; NAME:
;   cpbackup
;
; PURPOSE:
;   Copy a file to a backup file name.
;
; CALLING SEQUENCE:
;   cpbackup, filename
;
; INPUTS:
;   filename   - File to copy if it exists.
;
; OUTPUTS:
;
; COMMENTS:
;   Make a backup copy of the specified file by appending ".1", ".2", etc.
;   The first unused number is used as an appendix.
;
; EXAMPLES:
;
; BUGS:
;   This is only written to work with a Unix file system, since it spawns
;   the Unix "cp" command.
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   11-Mar-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

pro cpbackup, filename

   ;----------
   ; If the file FILENAME does not exist, then return

   fname1 = (findfile(filename, count=ct))[0]
   if (ct NE 1) then return

   ;----------
   ; Find the first backup file name that does not already exist.

   num = 0
   ct = 1
   while (ct NE 0) do begin
      num = num + 1
      backname = fname1 + '.' + strtrim(string(num),2)
      junk = (findfile(backname, count=ct))[0]
   endwhile

   ;----------
   ; Copy the file

   spawn, 'cp ' + fname1 + ' ' + backname

   return
end
