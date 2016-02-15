        FUNCTION PRODUCT, ARRAY, NAN = NAN
;+
; NAME:
;       PRODUCT()
; PURPOSE:
;       Calculate the product of all the elements of an array
; EXPLANATION:
;       This routine serves as an equivalent to the intrinsic PRODUCT() function
;       introduced in V5.6, except that the CUMULATIVE keyword is not available.
;
;       PRODUCT() is the multiplicative equivalent of TOTAL().   
; CALLING SEQUENCE:
;       Result = PRODUCT(ARRAY, [/NaN] )
; INPUT PARAMETERS:
;       ARRAY   = Array of elements to multiply together.  For instance, ARRAY
;                 could contain the dimensions of another array--then
;                 PRODUCT(ARRAY) would be the total number of elements of that
;                 other array.
; OUTPUT:
;       The result of the function is the total product of all the elements of
;       ARRAY.   If the input is double precision or 64bit integer, then the 
;       result will be the same; otherwise the result will be floating point.
; OPTIONAL OUTPUT KEYWORD:
;       /NAN - If set, then PRODUCT() will check for the presence of IEEE
;               floating point NaN values in the input array.
; RESTRICTIONS:
;       ARRAY must be a numerical type.

; PROCEDURE:
;      Vector multiplication in groups of powers of two make this operation
;      faster than a simple FOR loop.  The number of actual multiplications is 
;      still N_ELEMENTS(ARRAY).  Double precision should be used for the highest
;      accuracy when multiplying many numbers.
; MODIFICATION HISTORY:
;       William Thompson, Feb. 1992.
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Use vector algorithm from C. Markwardt's CMPRODUCT W. Landsman Nov. 2001   
;       Added /NAN keyword, documentation about V5.6 emulation W.L Nov. 2002       
;-
;
        ON_ERROR,2
;
;  Check the number of parameters.
;
        IF N_PARAMS() NE 1 THEN MESSAGE,'Syntax:  Result = PRODUCT(ARRAY)'
;
;  Check the type of ARRAY.
;
       TYPE = SIZE(ARRAY,/TYPE)
       IF TYPE EQ 0 THEN MESSAGE,'ARRAY not defined'
       IF TYPE EQ 7 THEN MESSAGE,'Operation illegal with string arrays'
       IF TYPE EQ 8 THEN MESSAGE,'Operation illegal with structures'
;
;  Calculate the product.    Make sure output is at least Floating pt.
;
       if (type EQ 5) or (type GE 14) then X = ARRAY ELSE X = FLOAT(ARRAY)
       if keyword_set(NAN) THEN BEGIN
               good = where( FINITE(X), NG)
               if NG GT 0 THEN X = X(WHERE(FINITE(X))) $
                          ELSE RETURN, !VALUES.F_NAN
       ENDIF
       N = N_ELEMENTS(X)
       WHILE N GT 1 DO BEGIN
           IF (N MOD 2) THEN X[0] = X[0] * X[N-1]
            N2 = FLOOR(N/2)
            X = X[0:N2-1] * X[N2:*]
            N = N2
        ENDWHILE
;
        RETURN, X[0]
        END
