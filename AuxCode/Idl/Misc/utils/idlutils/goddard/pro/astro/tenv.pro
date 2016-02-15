      FUNCTION tenv,dd,mm,ss
;+
; NAME:
;	TENV()
; PURPOSE:
;	Converts sexigesimal number or vector to decimal.  
; EXPLANATION:
;	Like TEN() but allows vector input.
;
; CALLING SEQUENCES:
;	Result = TENV( dd, mm )           ; result = dd + mm/60.
;	Result = TENV( dd, mm, ss)        ; result = dd + mm/60. + ss/3600.
;
; INPUTS:
;	dd - Sexigesimal element(s) corresponding to hours or degrees
;	mm - Sexigesimal element(s) corresponding to minutes
;	ss - Sexigesimal element(s) corresponding to seconds (optional)
;		The input parameters can be scalars or vectors.   However, the
;		number of elements in each parameter must be the same.
;
; OUTPUTS:
;	Result -  double, decimal equivalent of input sexigesimal 
;		quantities.  Same number of elements as the input parameters.
;		If the nth element in any of the input parameters is negative 
;		then the nth element in Result will also be negative.
;
; EXAMPLE:
;	If dd = [60,60,0], and mm = [30,-30,-30], then
;
;	IDL> Result = TENV(dd,mm)  ====>   Result =  [60.5,-60.5,-0.5]
;
; PROCEDURE:
;	Mostly involves checking arguments and setting the sign.
;
;   MODIFICATION HISTORY:
;	Written by W.B. Landsman           April, 1991
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_error,2                                 ;Return to caller

 npts = N_elements(dd)
 npar = N_params()
 if npts EQ 0 then begin
     print,'Syntax -  RESULT = TENV( dd, mm, ss)'
     return, 0.0d
 endif

 if ( npar EQ 1 ) then return,double( dd )   ;No need to check for neg values.

 value = double( abs(dd) ) 

 if ( npar GT 1 ) then begin               ;Add minutes/60., check for <0

      if N_elements(mm) NE npts then $
           message,'ERROR - Number of elements in each parameter must be equal'
      neg = (dd LT 0) OR  (mm LT 0) 
      value = value + abs(mm)/60.0d

 endif

 if ( npar GT 2 ) then begin               ;Add sec/3600., check for <0

      if N_elements(ss) NE npts then $
           message,'ERROR - Number of elements in each parameter must be equal'
      neg = neg OR  (ss LT 0) 
      value = value + abs(ss)/3600.0d

 endif

 neg = where( neg, Nfound )                  ;Account for negative values
 if ( Nfound GT 0 ) then value[neg] = -value[neg]

 return,value      
 end
