pro hprint, h, firstline
;+
; NAME:
;       HPRINT
; PURPOSE:
;       Print a FITS header (or other string array) one line at a time
; EXPLANATION:
;       The string array is  printed 1 line at a time.
;       Needed because IDL will add an extra space to the 80 character
;       FITS lines on TTY terminals, causing a space to appear between lines.
;
; CALLING SEQUENCE:
;       HPRINT, h, [ firstline ]
;
; INPUTS:
;       H - FITS header (or any other string array).
;
; OPTIONAL INPUT:
;       FIRSTLINE - scalar integer specifying the first line to begin 
;               displaying.   The default is FIRSTLINE = 1, i.e. display 
;               all the lines.     If Firstline is negative, then the first
;               line to be printed is counted backward from the last line.
;
; NOTES:
;       HPRINT has the following differences from the intrinsic PRINT procedure
;
;       (1) Arrays are printed one line at a time to avoid a space between 80
;               character lines
;       (2) Lines are trimmed with STRTRIM before being printed to speed up 
;               display
;       (3) The /more option is used for output. 
;
; EXAMPLE:
;       Read the header from a FITS file named 'test.fits' and display it at the
;       terminal beginning with line 50
;
;       IDL> h = headfits( 'test.fits')         ;Read FITS header
;       IDL> hprint, h, 50                      ;Display starting at line 50
;
;       To print the last 25 lines of the header
;
;       IDL> hprint, h, -25
;
; REVISION HISTORY:
;       Written W. Landsman                     July, 1990
;       Added test for user quit                July, 1991
;       Added optional FIRSTLINE parameter      November, 1992
;       Modified for when STDOUT is not a TTY W. Landsman  September 1995
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Fixed printing in IDLDE, C. Gehman      August, 1998
;       Skip PRINTF if IDL in demo mode  W. Landsman  October 2004
;       Fixed bug on non-terminals, William Thompson, 18-Oct-2004
;       
;-
  On_error,2                        ;Return to Caller
  compile_opt idl2

  if N_params() EQ 0 then begin
       print,'Syntax - HPRINT, h, [ firstline ]'
       return
  endif

  n = N_elements(h)
  if ( n EQ 0 ) then    $               ;Make sure input array is defined
     message,'String array (first parameter) not defined'

  if N_elements( firstline ) EQ 0 then firstline = 1
  if ( firstline[0] LT 0 ) then firstline = ( n + firstline[0]) > 1 < n  $
                                else firstline = firstline[0] > 1 < n

  if lmgr(/demo) then begin      ;in demo mode?
      for i=firstline-1, n-1 do print,h[i]
      return
  endif


; Now print the array one line at a time
  stdout = fstat(-1)
  if (stdout.isatty) and (not stdout.isagui) then begin  ;Open with /MORE if a TTY
      
      openw, outunit, filepath(/TERMINAL), /MORE, /GET_LUN  
      for i = firstline-1, n-1 do begin

     printf, outunit, strtrim( h[i] )
     if !ERR EQ 1 then goto, DONE      ;User entered "Q" in response to /more

  endfor

DONE:  free_lun, outunit

 endif else printf,-1,strtrim(h[firstline-1:*]), FORMAT='(A)'

  return
  end
