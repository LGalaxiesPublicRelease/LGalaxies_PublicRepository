;+
; NAME:
;   count_freelun
;
; PURPOSE:
;   Count the number of logical file pointers (LUNs) that are available
;
; CALLING SEQUENCE:
;   num = count_freelun()
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   num        - Number of logical file pointers (LUNs) that can be allocated
;                with GET_LUN before triggering an error
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   02-Mar-2004  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
function count_freelun

   ; As of IDL 6.0, one is still only able to allocate up to 29 file
   ; pointers total, numbered 100 to 128
   npossible = 30

   num = 0L
   lunarr = lonarr(npossible)
   catch, error_status
   if (error_status NE 0) then begin
      catch, /cancel
      !error_state.code = 0
      for j=0, num-1 do free_lun, lunarr[j]
      return, num
   endif

   for i=0L, npossible-1 do begin
      if (keyword_set(error_status) EQ 0) then get_lun, thislun
      lunarr[i] = thislun
      num = num + 1
   endfor

   for j=0, num-1 do free_lun, lunarr[j]
   return, num
end
;------------------------------------------------------------------------------
