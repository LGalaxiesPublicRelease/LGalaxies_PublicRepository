;+
; NAME:
;   gausspix1d
;
; PURPOSE:
;   Routine to evaluate a 1D gaussian function (plus polynomial terms)
;   integrated over a pixel of width 1.
;
; CALLING SEQUENCE:
;   gausspix, x, a, f, pder
;
; INPUTS:
;   x          - Dependent variable [NX].
;   a          - Coefficients:
;                A[0] = center of the Gaussian in X
;                A[1] = sigma width of the Gaussian
;                A[2] = normalization of the Gaussian; total area = A[1] * A[2]
;                A[3...] = polynomial coefficients for background terms
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   f          - Evaluated function [NX].
;   pder       - Array of partial derivatives [NX,NCOEFF].
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine can be passed as a function to the IDL CURVEFIT function.
;
;   If X is type DOUBLE, then the return values are also type DOUBLE;
;   otherwise they are type FLOAT.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   04-Aug-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro gausspix1d, x, a, f, pder

   ncoeff = n_elements(a)
   nx = n_elements(x)
   if (ncoeff LT 3) then $
    message, 'Too few elements for A'
   sqrtpi = sqrt(2 * !dpi)

   xlo = x - 0.5d
   xhi = x + 0.5d
   zlo = (xlo - a[0]) / a[1]
   zhi = (xhi - a[0]) / a[1]
   fint = a[2] * a[1] * (gaussint(zhi) - gaussint(zlo))
   f = fint

   if (ncoeff GT 3) then begin
      pterm = dblarr(nx,ncoeff)
      for icoeff=0, ncoeff-4 do begin
         pterm[*,icoeff] = (xhi^(icoeff+1) - xlo^(icoeff+1)) / (icoeff+1)
         f = f + a[icoeff+3] * pterm[*,icoeff]
      endfor
   endif

   if (n_params() GE 4) then begin
      pder = dblarr(nx, ncoeff)
      expzlo = exp(-zlo^2/2)
      expzhi = exp(-zhi^2/2)
      pder[*,0] = -a[2] * (expzhi - expzlo) / sqrtpi
      pder[*,1] = fint / a[1] - a[2] * (zhi * expzhi - zlo * expzlo) / sqrtpi
      pder[*,2] = fint / a[2]
      if (ncoeff GT 3) then $
       pder[*,3:ncoeff-1] = pterm
   endif

   if (size(x,/tname) NE 'DOUBLE') then begin
      f = float(f)
      if (n_params() GE 4) then pder = float(pder)
   endif

   return
end
;------------------------------------------------------------------------------
