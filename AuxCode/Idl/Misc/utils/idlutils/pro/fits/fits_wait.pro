;-----------------------------------------------------------------------
;+
; NAME:
;   fits_wait
;
; PURPOSE:
;   Wait for a FITS file to be fully written to disk.
;
; CALLING SEQUENCE:
;   qdone = fits_wait(filename, [deltat=, tmax=, _EXTRA=KeywordsForFitsRead ])
;
; INPUTS:
;   filename  - FITS file name.
;
; OPTIONAL INPUTS:
;   deltat    - Time to wait between attempts to check file; default to 10 sec.
;   tmax      - Maximum time to check file; default to 300 sec.
;   _EXTRA    - Keywords to pass to FITS_READ, such as /HEADER_ONLY or EXTEN=.
;
; OUTPUTS:
;   qdone     - Return 0 if the file was never read as a valid FITS file;
;               return 1 if it was.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The purpose of this routine is to be able to test when a FITS file
;   has been fully written to disk.  This is useful for real-time processes
;   that must check this.  The Goddard routine FITS_READ is used to
;   determine when a file is fully written.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   fits_read
;
; REVISION HISTORY:
;   02-Oct-2000  Written by David Schlegel, Princeton.
;-
;-----------------------------------------------------------------------
function fits_wait, filename, deltat=deltat, tmax=tmax, $
 _EXTRA=KeywordsForFitsRead

   ; Need at least 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - qdone = fits_wait(filename, [deltat=, tmax=, _EXTRA= ])'
      return, 0
   endif

   if (n_elements(deltat) EQ 0) then deltat = 10.0
   if (n_elements(tmax) EQ 0) then tmax = 300.0

   ;----------
   ; If there is a ".gz" extension, then assume that the file is OK
   ; and return 1, since I don't know how to test the validity of
   ; a gzipped FITS file w/out crashing.

   if (strmid(filename,(strlen(filename)-3)>0,3) EQ '.gz') then return, 1

   qdone = 0 ; =0 for not done, =1 for done

   tsum = 0.0
   while (qdone EQ 0 AND tsum LT tmax) do begin
      message = 0
      fits_read, filename, junk, /no_abort, message=message, $
       _EXTRA=KeywordsForFitsRead
      if (keyword_set(message)) then begin
         wait, deltat
         tsum = tsum + deltat
      endif else begin
         qdone = 1
      endelse
   endwhile

   return, qdone
end
;-----------------------------------------------------------------------
