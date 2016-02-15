
; Do a linear least squares fit of a 2-d parabola (6 free parameters)
; to an image with weights.
; model = a0 + a1 *x + a2*y + a3*x^2 + a4 * x * y + a5 * y^2
;

function fit_para2d, x, y, im, wt, model=model, chi2=chi2

     xbasis = fpoly(x,3)
     ybasis = fpoly(y,3)
  
     action = (x[*]*0) # replicate(1,6)
     action[*,0]  =1
     action[*,1] = xbasis[*,1]
     action[*,2] = ybasis[*,1]
     action[*,3] = xbasis[*,2]
     action[*,4] = xbasis[*,1] * ybasis[*,1]
     action[*,5] = ybasis[*,2]

     if (size(im))[0] EQ 1 then n_im = 1 $
     else n_im = (size(im))[2]
     res = fltarr(6, n_im)
     model = im*0.

     for i=0L, n_im-1 do begin
     
       a = action * (sqrt(wt[*,i]) # replicate(1,6))
       alpha = transpose(a) # a
       beta = transpose(action) # (im[*,i] * wt[*,i])

       if !version.release GE '5.6' then begin
         la_choldc, alpha, /double, status=status
         if status EQ 0 then  $
         res[*,i] = la_cholsol(alpha, beta, /double)
       endif else begin
         choldc, alpha, p, /double
         res[*,i] = cholsol(alpha, p, beta, /double)
       endelse
     endfor
     model = action # res
     chi2 = (im - model)^2 * wt  
     return, res
end

