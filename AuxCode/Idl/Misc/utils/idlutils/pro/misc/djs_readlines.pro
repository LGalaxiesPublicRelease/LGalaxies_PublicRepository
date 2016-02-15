;-----------------------------------------------------------------------
;+
; NAME:
;   djs_readlines()
;
; PURPOSE:
;   Read an ASCII file as one character string for each line.
;   If NHEAD is specified and greater than zero, then that number
;   of lines is read in first as header strings in HEAD.
;
; CALLING SEQUENCE:
;   Data = djs_readlines( infile, [ nhead=nhead, Head=Head ] )
;
; INPUT:
;   infile:      Input file name
;
; OPTIONAL INPUT:
;   nhead:       Number of lines in header
;
; OUTPUTS:
;   Data:        Array of strings, one for each data line
;
; OPTIONAL OUTPUTS:
;   Head:        Array of strings, one for each header line
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Written by D. Schlegel, 29 May 1997, Durham
;-
;-----------------------------------------------------------------------
function djs_readlines, infile, nhead=nhead, Head=Head

   ; Need 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - Data = djs_readlines( infile, [ nhead=nhead, Head=Head ] )'
      return, -1
   endif

   if (NOT keyword_set(nhead)) then nhead = 0
   nline = numlines(infile)
   ndata = nline - nhead
   if (nline LT nhead) then return, -1

   ; Open file
   get_lun, ilun
   openr, ilun, infile
   tmpString = ''

   ; Read header lines
   if (nhead GT 0) then begin
      Head = strarr(nhead)
      for iline=0L, nhead-1 do begin
         readf, ilun, tmpString
         Head[iline] = tmpString
      endfor
   endif

   ; Read data lines
   if (ndata LE 0) then return, -1
   Data = strarr(ndata)
   for iline=0L, ndata-1 do begin
      readf, ilun, tmpString
      Data[iline] = tmpString
   endfor

   ; Close data file
   close, ilun
   free_lun, ilun

   return, Data
end 
;-----------------------------------------------------------------------
