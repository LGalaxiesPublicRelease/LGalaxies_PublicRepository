      FUNCTION sixty,scalar
;+
; NAME:
;	SIXTY()
; PURPOSE:
;	Converts a decimal number to sexigesimal.
; EXPLANATION:
;	Reverse of the TEN() function.
;
; CALLING SEQUENCE:
;	X = SIXTY( SCALAR ) 
;
; INPUTS:
;	SCALAR -- Decimal quantity.  
; OUTPUTS:
;	Function value returned = floating real vector of three elements, 
;	sexigesimal equivalent of input decimal quantity.
;	A negative number is signified by making the first non-zero
;	element of the output vection negative.
;
; PROCEDURE:
;	Mostly involves checking arguments and setting the sign.
;
; EXAMPLE:
;	If x = -0.345d then sixty(x) = [0.0, -20.0, 42.0]
;
; MODIFICATION HISTORY:
;	Written by R. S. Hill, STX, 19-OCT-87         
;	Output changed to single precision.  RSH, STX, 1/26/88
;	Accept single element vector   W. Landsman   Sep. 1996
;	Converted to IDL V5.0   W. Landsman   September 1997
;-

      if N_elements(scalar) NE 1 then begin
	      message,'ERROR - First parameter must contain 1 element',/CON
	      return,replicate(100.0e0,3)
      endif	

      ss=abs(3600.0d0*scalar)
      mm=abs(60.0d0*scalar) 
      dd=abs(scalar) 
      result=fltarr(3)
      result[0]=float(fix(dd))
      result[1]=float(fix(mm-60.0d0*result[0]))
      result[2]=float(ss-3600.d0*result[0]-60.0d0*result[1])
      if scalar[0] lt 0.0d0 then begin 
         if result[0] ne 0 then result[0] = -result[0] else $
         if result[1] ne 0 then result[1] = -result[1] else $
         result[2] = -result[2]
      endif

      return,result
      end
