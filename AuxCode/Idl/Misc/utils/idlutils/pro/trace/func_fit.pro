;+
; NAME:
;   func_fit
;
; PURPOSE:
;   Fit X, Y positions to a functional form.
;
; CALLING SEQUENCE:
;   res = func_fit( x, y, ncoeff, [invvar=, function_name=, $
;    ia=, inputans=, yfit=, double=double, _EXTRA= ]
;
; INPUTS:
;   x          - X values (independent variable)
;   y          - Y values (dependent variable)
;   ncoeff     - Number of coefficients to fit
;   invvar     - Weight values (inverse variance)
;
; OPTIONAL KEYWORDS:
;   function_name - Function to fit; options are:
;                'legendre'
;                'chebyshev'
;                Default to 'legendre'
;   ia         - Array specifying free (1) and fixed (0) variables [NCOEFF]
;   inputans   - If holding parameters fixed, set this array to those values
;                [NCOEFF]
;   double     - If set, or if X, Y, or INVVAR are double-precision, then
;                return double-precision values
;   _EXTRA     - Optional keywords for fitting function
;
; OUTPUTS:
;   res        - Fit coefficients [NCOEFF]
;
; OPTIONAL OUTPUTS:
;   yfit       - Fit evaluated at the points X
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   fchebyshev()
;   flegendre()
;   fpoly()
;
; REVISION HISTORY:
;   10-Sep-1999  Written by S. Burles
;   16-Nov-1999  Modified by D. Schlegel to never fit more coefficients
;                than there are data points.
;   10-Jul-2001  Added fpoly, S. Burles
;-
;------------------------------------------------------------------------------
function func_fit, x, y, ncoeff, invvar=invvar, function_name=function_name, $
 ia=ia, inputans=inputans, yfit=yfit, inputfunc=inputfunc, $
 double=double1, _EXTRA=KeywordsForFunc

   if (N_params() LT 3) then begin
      print,'res = func_fit( x, y, ncoeff, [ invvar=, function_name='
      print,' ia=, inputans=, yfit=, inputfunc=, /double ] )'
      return, 0
   endif

   if (size(x,/tname) EQ 'DOUBLE' OR size(y,/tname) EQ 'DOUBLE' $
    OR size(invvar,/tname) EQ 'DOUBLE' OR keyword_set(double1)) then $
    double = 1B

   if (NOT keyword_set(function_name)) then function_name = 'flegendre'
   nx = N_elements(x)
   ny = N_elements(y)
   if (nx NE ny) then message, 'Dimensions of X and Y do not agree'

   if (keyword_set(invvar)) then begin
      nw = N_elements(invvar)
      if (nx NE nw) then message, 'Dimensions of X and INVVAR do not agree'
   endif else begin
      invvar = intarr(nx) + 1
   endelse

   ; If IA is not set, then let all coefficients be fit (e.g., IA[*]=1)
   if (n_elements(ia) NE ncoeff) then ia = bytarr(ncoeff) + 1

   goodia = ia NE 0 ; Indices of coefficients to fit

   ; If coefficients are being fixed but INPUTANS is not specified,
   ; then declare INPUTANS and set all those fixed coefficients equal to zero.
   fixed = where(goodia EQ 0, nfix)
   if (nfix GT 0) then $
    if (NOT keyword_set(inputans)) then $
     inputans = keyword_set(double) ? dblarr(ncoeff) : fltarr(ncoeff)

   ;------
   ; Select unmasked points

   igood = where(invvar GT 0, ngood)
   res = keyword_set(double) ? dblarr(ncoeff) : fltarr(ncoeff)

   if (ngood EQ 0) then begin

      ; Case with no good data points.
      ; Set all fit coefficients equal to zero.
      ; Also, set all the YFIT values equal to zero.

      yfit = 0 * y

   endif else if (ngood EQ 1) then begin

      ; Case with only 1 good data point.
      ; Set the first fit coefficient equal to the value of that one data point
      ; (and other coefficients equal to zero).
      ; Also, set all the YFIT values to that same value.

      res[0] = y[igood[0]]
      yfit = 0 * y
      yfit[*] = y[igood[0]]

   endif else begin

      ; Do not fit more coefficients than there are data points
      ncfit = ncoeff < ngood

      if (function_name EQ 'poly') then legarr = fpoly(x, ncfit)
      if (function_name EQ 'legendre' OR function_name EQ 'flegendre') then $
       legarr = flegendre(x, ncfit)
      if (function_name EQ 'chebyshev' OR function_name EQ 'fchebyshev') then $
       legarr = fchebyshev(x, ncfit)
;      if (function_name EQ 'flegendre') then $
;       legarr = flegendre(x, ncfit)
;      if (function_name EQ 'fchebyshev') then $
;       legarr = fchebyshev(x, ncfit)
      if (function_name EQ 'fchebyshev_split') then $
       legarr = fchebyshev_split(x, ncfit)

      if keyword_set(inputfunc) then $
        legarr = legarr * (inputfunc # replicate(1,ncfit))

      ; Subtract fixed terms first
     
      yfix = y * 0.0 
      nonfix = where(goodia[0:ncfit-1], nparams)
      fixed = where(goodia[0:ncfit-1] EQ 0, nfix)
      
      if (nfix GT 0) then begin
         yfix = legarr # (inputans * (1 - goodia)) 
         ysub = y - yfix
         finalarr = legarr[*, nonfix]
      endif else begin
         finalarr = legarr
         ysub = y
      endelse

      tempvec = keyword_set(double) ? dblarr(nparams) : fltarr(nparams)
      extra2 = finalarr * ((tempvec + 1) ## (invvar > 0))
      alpha = transpose(finalarr) # extra2

      ; Don't send just one parameter to SVD fit

      if (nparams GT 1) then begin
         beta = transpose( (ysub * (invvar > 0)) # finalarr)
         svdc, alpha, w, u, v, /double
         res[nonfix] = svsol(u, w, v, beta, /double)
      endif else begin
         res[0] = total(ysub * invvar * finalarr) / alpha
      endelse

      if (nfix GT 0) then res[fixed] = inputans[fixed]

      yfit = legarr # res[0:ncfit-1]

   endelse

   return, res
end
