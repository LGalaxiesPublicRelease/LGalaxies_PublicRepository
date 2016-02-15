;------------------------------------------------------------------------------
;+
; NAME:
;   djs_hex2bin
;
; PURPOSE:
;   Convert hexadecimal number(s) to binary numbers.
;
; CALLING SEQUENCE:
;   binval = djs_hex2bin(hexval, [ndigit=ndigit])
;
; INPUTS:
;   hexval:   String or array of strings containing hexadecimal number(s)
;
; OPTIONAL INPUTS:
;   ndigit:   Number of binary digits in output; if not supplied, then the
;             minimum number of digits are used
;
; OUTPUTS:
;   intval:   Byte array(s) of binary values, dimensioned intval(ndigit,nnum)
;
; EXAMPLE:
;   PRINT, DJS_HEX_TO_BINARY( '1a' )
;   IDL prints the result
;     0   1   0   1   1
;
;   One can truncate to only the 4 least significan digits by setting NDIGIT:
;   PRINT, DJS_HEX_TO_BINARY( '1a', NDIGIT=4 )
;   IDL prints the result
;     0   1   0   1
;
;   PRINT, DJS_HEX_TO_BINARY( ['1a', '2b', '3'] )
;   IDL prints the result
;     0   1   0   1   1   0
;     1   1   0   1   0   1
;     1   1   0   0   0   0
;
; PROCEDURES CALLED:
;   djs_hex2int()
;   djs_int2bin()
;
; REVISION HISTORY:
;   Written D. Schlegel, 30 June 1997, Durham
;-
;------------------------------------------------------------------------------
function djs_hex2bin, hexval, ndigit=ndigit
 
   ; Need 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - binval = djs_hex2bin( hexval, [ndigit=ndigit] )'
      return, -1
   endif

   intval = djs_hex2int(hexval)
   binval = djs_int2bin(intval, ndigit=ndigit)

   return, binval
end
;------------------------------------------------------------------------------
