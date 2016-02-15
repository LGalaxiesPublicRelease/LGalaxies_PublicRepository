;------------------------------------------------------------------------------
;+
; NAME:
;   djs_hex2int
;
; PURPOSE:
;   Convert hexadecimal number(s) to long integer(s).
;
; CALLING SEQUENCE:
;   intval = djs_hex2int(hexval)
;
; INPUTS:
;   hexval:   String or array of strings containing hexadecimal number(s)
;
; OUTPUTS:
;   intval:   Long integer or array or long integers
;
; EXAMPLE:
;   PRINT, DJS_HEX_TO_INT( '1a' )
;   IDL prints the result
;     26
;
;   PRINT, DJS_HEX_TO_INT( ['1a', '2b', '3'] )
;   IDL prints the result
;     26  43   3
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Written D. Schlegel, 30 June 1997, Durham
;-
;------------------------------------------------------------------------------
function djs_hex2int, hexval
 
   ; Need 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - intval = djs_hex2int( hexval )'
      return, -1
   endif

   ndigit = strlen(hexval)
   maxdigit = max(ndigit)
   mfac = 16L^indgen(maxdigit)

   ; Case where input value is one string
   if ((size(hexval))(0) EQ 0) then begin
      bval = byte(hexval)
      cval = (bval LT 58) * (bval-48) + (bval GT 96) * (bval-87)
      intval = long( total( reverse(cval) * mfac) )

   ; Case where input value is an array of strings
   endif else begin
      nnum = (size(hexval))(1)
      intval = lonarr(nnum)
      for inum=0L, nnum-1 do begin
         bval = byte(hexval(inum))
         cval = (bval GE 48 AND bval LE 57) * (bval-48) $
              + (bval GE 65 AND bval LE 70) * (bval-55) $
              + (bval GE 97 AND bval LE 102) * (bval-87)
         intval(inum) = long( total( reverse(cval) * mfac) )
      endfor
   endelse

   return, intval
end
;------------------------------------------------------------------------------
