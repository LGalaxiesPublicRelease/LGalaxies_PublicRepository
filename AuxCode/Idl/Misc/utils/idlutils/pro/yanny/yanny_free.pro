;+
; NAME:
;   yanny_free
;
; PURPOSE:
;   Free memory allocated from reading a Yanny parameter file
;
; CALLING SEQUENCE:
;   yanny_free, pdata
;
; INPUTS:
;   pdata      - Array of pointers to all strucutures read by YANNY_READ.
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
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
;   15-Nov-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro yanny_free, pdata

   if (N_params() LT 1) then begin
      print, 'Syntax - yanny_free, pdata'
      return
   endif

   if (NOT keyword_set(pdata)) then return

   for i=0, N_elements(pdata)-1 do begin
      ptr_free, pdata[i]
   endfor

   pdata = 0

   return
end
;------------------------------------------------------------------------------
