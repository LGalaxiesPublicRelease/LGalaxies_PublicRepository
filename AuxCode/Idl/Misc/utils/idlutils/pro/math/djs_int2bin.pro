;------------------------------------------------------------------------------
;+
; NAME:
;   djs_int2bin
;
; PURPOSE:
;   Convert integer number(s) to binary numbers.
;
; CALLING SEQUENCE:
;   binval = djs_int2bin(val, [ndigit=ndigit])
;
; INPUTS:
;   val:      Integer number(s)
;
; OPTIONAL INPUTS:
;   ndigit:   Number of binary digits in output; if not supplied, then the
;             minimum number of digits are used
;
; OUTPUTS:
;   binval:   Byte array(s) of binary values
;
; PROCEDURES CALLED:
;   djs_floor()
;
; REVISION HISTORY:
;   Written D. Schlegel, 30 June 1997, Durham
;   31-Jul-1998  DJS - Subscripts modified to IDL 5 convention.
;   03-Aug-1999  DJS - Modified to work with signed integers by
;                first converting to unsigned integers.
;-
;------------------------------------------------------------------------------
function djs_int2bin, val, ndigit=ndigit
 
   ; Need 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - binval = djs_int2bin( val, [ndigit=ndigit] )'
      return, -1
   endif

   ; If there is a negative value, then change to an unsigned integer
   tname = size(val,/tname)
   case tname of
   'BYTE': begin
      intval = byte(val)
      TWO = byte(2)
      end
   'INT': begin
      intval = uint(val)
      TWO = uint(2)
      end
   'UNIT': begin
      intval = uint(val)
      TWO = uint(2)
      end
   'LONG': begin
      intval = ulong(val)
      TWO = ulong(2)
      end
   'ULONG': begin
      intval = ulong(val)
      TWO = ulong(2)
      end
   'LONG64': begin
      intval = ulong64(val)
      TWO = ulong64(2)
      end
   'ULONG64': begin
      intval = ulong64(val)
      TWO = ulong64(2)
      end
   else: begin
      message, 'Unsupported type for VAL'
      end
   endcase

   ; The addition of 0.01 below is just to prevent rounding-down problems
   if (max(intval) LE 0) then maxdigit = 1 $
   else maxdigit = djs_floor( alog(max(intval)) / alog(2.0) + 0.01 ) + 1

   mfac = TWO^indgen(maxdigit)

   ; Case where input value is a scalar
   if ((size(intval))[0] EQ 0) then begin
      nnum = 1
      binval = bytarr(maxdigit)
      for idig = 0, maxdigit-1 do $
       binval[idig] = (intval AND TWO^idig) NE 0

   ; Case where input value is a vector
   endif else begin
      nnum = (size(intval))[1]
      binval = bytarr(maxdigit, nnum)
      for idig = 0, maxdigit-1 do $
       binval[idig,*] = (intval AND TWO^idig) NE 0
   endelse

   if (keyword_set(ndigit)) then begin
      if (ndigit LT maxdigit) then $
       binval = binval[0:ndigit-1,*]
      if (ndigit GT maxdigit) then $
       binval = [binval, bytarr(ndigit-maxdigit,nnum)]
   endif

   return, binval
end
;------------------------------------------------------------------------------
