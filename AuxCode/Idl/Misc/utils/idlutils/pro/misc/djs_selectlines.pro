;-----------------------------------------------------------------------
;+
; NAME:
;   djs_selectlines
;
; PURPOSE:
;   Select the line numbers specified by INDX of a file, and print either
;   to the standard output or to another file.
;
;   This is not yet optimized for very large files, as it will read in
;   all of the requested lines (though not all the lines) into memory first.
;
; CALLING SEQUENCE:
;   djs_selectlines, infile, [ indx=indx, nhead=nhead, outfile=outfile ]
;
; INPUTS:
;   infile:      Input file name
;
; OPTIONAL INPUTS:
;   indx:        Array of line numbers to select (0-indexed); default to all.
;                The indices start with 0 for the first data line.
;   nhead:       Number of lines in header
;   outfile:     Output file name; if not set then print to terminal
;
; OUTPUTS:
;
; PROCEDURES CALLED:
;   djs_readilines()
;
; REVISION HISTORY:
;   Written by D. Schlegel, 25 September 1997, Durham
;-
;-----------------------------------------------------------------------
pro djs_selectlines, infile, indx=indx, nhead=nhead, outfile=outfile

   ; Need at least 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - djs_selectlines, infile, [ indx=indx, nhead=nhead,'
      print, ' outfile=outfile ]'
      return
   endif

   ; Read selected lines
   Data = djs_readilines( infile, indx=indx, nhead=nhead, Head=Head )

   ; Open file
   if (keyword_set(outfile)) then begin
      get_lun, olun
      openw, olun, outfile
   endif else begin
      olun = -1
   endelse

   ; Print the header lines
   if (nhead GT 0) then begin
      printf, olun, Head, format='(a)'
   endif

   ; Print the data lines
   if (N_elements(Data) GT 0) then begin
      printf, olun, Data, format='(a)'
   endif

   ; Close data file
   if (keyword_set(outfile)) then begin
      close, olun
      free_lun, olun
   endif

   return
end 
;-----------------------------------------------------------------------
