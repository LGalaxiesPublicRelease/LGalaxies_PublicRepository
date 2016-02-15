;+
; NAME:
;   rmfile
;
; PURPOSE:
;   Delete file from disk.
;
; CALLING SEQUENCE:
;   rmfile, filename
;
; INPUTS:
;   filename   - File to delete.
;
; OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   14-Oct-1999  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------
pro rmfile, filename

   if (NOT keyword_set(filename)) then return

   ;----------
   ; Call this routine recursively if FILENAME is an array
   nfile = n_elements(filename)
   if (nfile GT 1) then begin
      for ifile=0, nfile-1 do rmfile, filename[ifile]
      return
   endif

   ; It can be the case that one element of the array is blank.
   if (NOT keyword_set(filename)) then return

   get_lun, ilun
   openr, ilun, filename, /delete, error=err
   if (err NE 0) then message, !err_string, /informational
   close, ilun
   free_lun, ilun

   return
end
