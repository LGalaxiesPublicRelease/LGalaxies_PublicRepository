
;+
; NAME:
;   djs_sfit
;
; PURPOSE:
;   Surface-fitting code to tabulated data (optionally with errors).
;
; CALLING SEQUENCE:
;   acoeff = djs_sfit( fval, xval, yval, degreex, degreey, $
;    [ sqivar=, yfit= ] )
;
; INPUTS:
;   fval       - Function values at XVAL,YVAL.
;   xval       - X coordinate values
;   yval       - Y coordinate values
;   degreex    - Degree of polynomial fit in X; 1 for linear, 2 for quadratic
;   degreey    - Degree of polynomial fit in Y; 1 for linear, 2 for quadratic
;
; OPTIONAL INPUTS:
;   sqivar     - Inverse sigma, which are the weights
;
; OUTPUTS:
;   acoeff     - Fit coefficients as [DEGREEX+1,DEGREEY+1] array
;
; OPTIONAL OUTPUTS:
;   yfit       - Fit values
;
; COMMENTS:
;
; EXAMPLES:
;   Create a random 2-dimensional field with a gradient in the X direction,
;   and fit to a quadratic function in both X and Y:
;     IDL> xval = dindgen(100) # replicate(1,100) / 100.
;     IDL> yval = replicate(1,100) # dindgen(100) / 100.
;     IDL> image = smooth(randomu(1234,100,100),11,/edge) + 0.2*xval^2
;     IDL> acoeff = djs_sfit(image,xval,yval,2,2,yfit=yfit)
;   Display the original image, and then the residual between that
;   image and the fit:
;     IDL> atv, image
;     IDL> atv, image - yfit
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINS:
;
; REVISION HISTORY:
;   25-Oct-2002  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function djs_sfit, fval, xval, yval, degreex, degreey, $
 sqivar=sqivar, yfit=yfit

   npts = n_elements(fval)
   if (n_params() LT 5) then $
    message, 'Wrong number of parameters'
   if (n_elements(xval) NE npts OR n_elements(yval) NE npts) then $
    message, 'Dimensions of FVAL,XVAL,YVAL must agree'
   if (degreex LT 0 OR degreey LT 0) then $
    message, 'DEGREEX,DEGREEY must be non-negative integers'

   ncoeff = (degreex+1) * (degreey+1)
   if (keyword_set(sqivar)) then thisvec = reform(double(sqivar),npts) $
    else thisvec = dblarr(npts)+1.d0
   if (n_elements(thisvec) NE npts) then $
    message, 'Number of elements in SQIVAR must agree with FVAL'
   thisx = reform(double(xval), npts)
   thisy = reform(double(yval), npts)

   mmatrix = dblarr(npts, ncoeff, /nozero)
   for xp=0, degreex do begin
      for yp=0, degreey do begin
         tmpvec = thisvec
         if (xp GT 0) then tmpvec = tmpvec * thisx^xp
         if (yp GT 0) then tmpvec = tmpvec * thisy^yp
         mmatrix[*,yp*(degreex+1)+xp] = tmpvec
      endfor
   endfor

   bvec = reform(double(fval),npts) * thisvec
   if (ncoeff EQ 1) then begin
      mt_m = total(mmatrix * mmatrix)
      mt_b = total(mmatrix * bvec)
      mmi = 1. / (mt_m + (mt_m EQ 0))
   endif else begin
      mt_m = matrix_multiply(mmatrix, mmatrix, /atranspose)
      mt_b = matrix_multiply(mmatrix, bvec, /atranspose)
      mmi = invert(mt_m, /double)
   endelse
   acoeff = mmi # mt_b

   if (arg_present(yfit)) then $
    yfit = reform(acoeff ## mmatrix, size(yval, /dimens))

   return, reform(acoeff, degreex+1, degreey+1)
end
;------------------------------------------------------------------------------
