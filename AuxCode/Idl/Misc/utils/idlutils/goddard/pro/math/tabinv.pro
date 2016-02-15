PRO TABINV, XARR, X, IEFF, FAST = fast
;+ 
; NAME:
;       TABINV     
; PURPOSE:  
;       To find the effective index of a function value in an ordered vector.
;       
; CALLING SEQUENCE:
;       TABINV, XARR, X, IEFF, [/FAST]
; INPUTS:
;       XARR - the vector array to be searched, must be monotonic
;               increasing or decreasing
;       X    - the function value(s) whose effective
;               index is sought (scalar or vector)
;
; OUTPUT:
;       IEFF - the effective index or indices of X in XARR
;               real or double precision, same # of elements as X
;
; OPTIONAL KEYWORD INPUT:
;       /FAST - If this keyword is set, then the input vector is not checked
;               for monotonicity, in order to improve the program speed.
; RESTRICTIONS:
;       TABINV will abort if XARR is not monotonic.  (Equality of 
;       neighboring values in XARR is allowed but results may not be
;       unique.)  This requirement may mean that input vectors with padded
;       zeroes could cause routine to abort.
;
; PROCEDURE:
;       VALUE_LOCATE() is used to find the values XARR[I]
;       and XARR[I+1] where XARR[I] < X < XARR[I+1].
;       IEFF is then computed using linear interpolation 
;       between I and I+1.
;               IEFF = I + (X-XARR[I]) / (XARR[I+1]-XARR[I])
;       Let N = number of elements in XARR
;               if x < XARR[0] then IEFF is set to 0
;               if x > XARR[N-1] then IEFF is set to N-1
;
; EXAMPLE:
;       Set all flux values of a spectrum (WAVE vs FLUX) to zero
;       for wavelengths less than 1150 Angstroms.
;         
;       IDL> tabinv, wave, 1150.0, I
;       IDL> flux[ 0:fix(I) ] = 0.                         
;
; FUNCTIONS CALLED:
;       None
; REVISION HISTORY:
;       Adapted from the IUE RDAF                     January, 1988         
;       More elegant code  W. Landsman                August, 1989
;       Mod to work on 2 element decreasing vector    August, 1992
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Updated for V5.3 to use VALUE_LOCATE()     W. Landsman January 2000
;       Work when both X and Xarr are integers     W. Landsman August 2001
;-               
 On_error,2

 if N_params() LT 3 then begin
     print,'Syntax- TABINV, XARR, X, I, [/FAST]'
     return
 endif

 Npoints = N_elements(xarr) & npt= npoints - 1

 if not keyword_set(fast) then begin

 if ( Npoints LE 1 ) then message, /TRACE, $
   'Search vector (first parameter) must contain at least 2 elements'

 ; Test for monotonicity 

 i = xarr - shift( xarr,1)
 i = i[1:*]               ;Added 15-Aug-1992 to properly interpret 2 element
 a = where( i GE 0, N)    ;decreasing vector

 if ( N NE npt) then begin 

     a = where(i LE 0, N)  ; Test for decreasing array
     if ( N NE npt ) then  message, /TRACE, $
       'ERROR - First parameter must be a monotonic vector' 
 endif
 endif

 
 ieff = float(VALUE_LOCATE(xarr,x)) 
 g = where( (ieff LT npt) and (ieff GE 0), Ngood)
 if Ngood GT 0 then begin
      neff = ieff[g]
     diff = x[g] - xarr[neff]
     if size(diff,/TNAME) NE 'DOU' then diff = float(diff) 
     ieff[g] = neff +  diff / (xarr[neff+1] - xarr[neff] ) 
 endif
 ieff = ieff > 0.0

 return
 end
