;+
; NAME:
;   poly_iter
;
; PURPOSE:
;   Calls IDL poly_fit iteratively with outlier rejection
;
; CALLING SEQUENCE:
;   poly_iter, x, y, ndeg, nsig, yfit, coeff=coeff
;
; INPUTS:
;   x, y    - indep, dep variables
;   ndeg    - degree of polynomial 
;   yfit    - fit y at given x values
;
; OUTPUTS:
;   coeff   - array of coefficients
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   poly_fit
;
; REVISION HISTORY:
;   moved from hoggpt 10-Jan-2002
;-
;------------------------------------------------------------------------------
PRO poly_iter, x, y, ndeg, nsig, yfit, coeff=coeff

  good = bytarr(n_elements(x))+1B
  w = lindgen(n_elements(x))

  FOR i=1, 5 DO BEGIN 
      if (!version.release LT '5.4') then $
       coeff = poly_fit(x[w], y[w], ndeg, yfit) $
      else $
       coeff = poly_fit(x[w], y[w], ndeg, yfit, /double)

      res = y[w]-yfit
      sig = stddev(res, /double)
      good[w] = good[w]*(abs(res) LE nsig*sig)
      w = where(good)
  ENDFOR 

  if (!version.release LT '5.4') then $
   coeff = poly_fit(x[w], y[w], ndeg) $
  else $
   coeff = poly_fit(x[w], y[w], ndeg, /double)


  yfit = poly(x, coeff)

  return
END
