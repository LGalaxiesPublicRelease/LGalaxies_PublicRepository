;-----------------------------------------------------------------------
;+
; NAME:
;   djs_readilines()
;
; PURPOSE:
;   Read selected lines of an ASCII file as one character string for each line.
;   If NHEAD is specified and greater than zero, then that number
;   of lines is read in first as header strings in HEAD.
;
; CALLING SEQUENCE:
;   Data = djs_readilines( infile, indx=indx, [ nhead=nhead, Head=Head ] )
;
; INPUT:
;   infile:      Input file name
;
; OPTIONAL INPUT:
;   nhead:       Number of lines in header
;   indx:        Line numbers to read (0-indexed); if not set, then
;                default to reading all data lines in their current order.
;                The indices start with 0 for the first data line.
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
function djs_readilines, infile, indx=indx, nhead=nhead, Head=Head

   ; Need 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - Data = djs_readilines( infile, indx=indx,'
      print, ' [ nhead=nhead, Head=Head ] )'
      return, -1
   endif

   if (NOT keyword_set(nhead)) then nhead = 0
   nline = numlines(infile)
   ndata = nline - nhead
   if (nline LT nhead) then return, -1
   if (n_elements(indx) EQ 0L) then indx = lindgen(ndata)

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

   ; Make a sorted list of indices
   sind = sort(indx)

   ; Construct output data array
   nout = N_elements(indx)
   Data = strarr(nout)

   ; Read data lines
   iline = 0
   if (ndata LE 0) then return, -1
   for iout=0L, nout-1 do begin
      if (iout EQ 0) then nread = indx[sind[0]] + 1 $
      else nread = indx[sind[iout]] - indx[sind[iout-1]]
      if (nread GT 0) then begin
         for iread=0L, nread-1 do $
          readf, ilun, tmpString
      endif
      Data[sind[iout]] = tmpString
   endfor

   ; Close data file
   close, ilun
   free_lun, ilun

   return, Data
end 
;-----------------------------------------------------------------------
