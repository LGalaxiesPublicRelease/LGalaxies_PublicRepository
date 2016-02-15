      FUNCTION ten,dd,mm,ss
;+
; NAME:
;	TEN()
; PURPOSE:
;	Converts a sexigesimal number to decimal.
; EXPLANATION:
;	Inverse of the SIXTY() function.
;
; CALLING SEQUENCES:
;	X = TEN( [ HOUR_OR_DEG, MIN, SEC ] )
;	X = TEN( HOUR_OR_DEG, MIN, SEC )
;	X = TEN( [ HOUR_OR_DEG, MIN ] )
;	X = TEN( HOUR_OR_DEG, MIN )
;	X = TEN( [ HOUR_OR_DEG ] )      <--  Trivial cases
;	X = TEN( HOUR_OR_DEG )        <--
;
; INPUTS:
;	HOUR_OR_DEG,MIN,SEC -- Scalars giving sexigesimal quantity in 
;		in order from largest to smallest.    
;
; OUTPUTS:
;	Function value returned = double real scalar, decimal equivalent of
;	input sexigesimal quantity.  A minus sign on any nonzero element
;	of the input vector causes all the elements to be taken as
;	< 0.
;
; PROCEDURE:
;	Mostly involves checking arguments and setting the sign.
;
;	The procedure TENV can be used when dealing with a vector of 
;	sexigesimal quantities.
;
; MODIFICATION HISTORY:
;	Written by R. S. Hill, STX, 21 April 87       
;	Modified to allow non-vector arguments.  RSH, STX, 19-OCT-87
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
      np = N_params()

      if (np eq 1) then begin            
         vector=dd
      endif else begin
         if (np lt 1) or (np gt 3) then goto,bad_args
         vector=dblarr(3)
         vector[0]=dd
         vector[1]=mm
         if np gt 2 then vector[2]=ss
      endelse
      sz = size(vector)
      ndim = sz[0]
      if (ndim eq 0) then return,double(vector)
      facs=[1.0d0,60.0d0,3600.0d0]
      nel = sz[1]
      sign = +1.0d0
      signflag = total(vector lt dblarr(nel))
      if (signflag gt 0.0d0) then sign = -1.0d0
      vector = abs(vector)
      decim = double(vector[0])
      i = 1
      while (i le nel-1) do begin
         decim = decim + double(vector[i])/facs[i]
         i = i + 1
      endwhile
      return,decim*sign
bad_args:    
      print,'Argument(s) should be hours/degrees, minutes (optional),'
      print,'seconds (optional)   in vector or as separate arguments.'
      print,'If any one number negative, all taken as negative.'
      return,0.0d0     
      end
