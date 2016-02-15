;
;  An internal routine for bspline_extract and bspline_radial
;   which solve a general banded correlation matrix which is represented
;   by the variable "action".  This routine only solves the lienar system
;   once, and stores the coefficients in sset.
;  A non-zero return value signifies a failed inversion
;
;
;------------------------------------------------------------------------------
function bspline_workit, xdata, ydata, invvar, action, sset, $
            lower=lower, upper=upper, npoly=npoly, nord=nord, yfit=yfit, $
            covariance=covariance, alpha=alpha

   nord = sset.nord
   npoly = sset.npoly
   goodbk = where(sset.bkmask[nord:*] NE 0, nbkpt)

   if (nbkpt LT nord) then begin
      if (arg_present(yfit)) then yfit = fltarr(n_elements(ydata))
      return, -2L
   endif

   nn = nbkpt 
   nfull = nn * npoly
   bw = npoly * nord   ; this is the bandwidth

   ;  The next line is REQUIRED to fill a1

   a2 = action * sqrt(invvar # replicate(1,bw))

   alpha = dblarr(bw,nfull+bw)
   beta = dblarr(nfull+bw)

   bi = lindgen(bw)
   bo = lindgen(bw)
   for i=1L, bw-1 do bi = [bi, lindgen(bw-i)+(bw+1)*i]
   for i=1L, bw-1 do bo = [bo, lindgen(bw-i)+bw*i]

   for i=0L, nn-nord do begin

      itop = i * npoly
      ibottom = (itop < nfull + bw) - 1
       
      ict = upper[i] - lower[i] + 1
  
      if (ict GT 0) then begin

         work = a2[lower[i]:upper[i],*] ## transpose(a2[lower[i]:upper[i],*])
         wb   =  (ydata[lower[i]:upper[i]]*sqrt(invvar[lower[i]:upper[i]])) $
                            # a2[lower[i]:upper[i],*] 

         alpha[bo+itop*bw] = alpha[bo+itop*bw] + work[bi]
         beta[itop:ibottom] = beta[itop:ibottom] + wb

      endif
   endfor

   ; Drop break points where minimal influence is located

   min_influence = 1.0e-10 * total(invvar) / nfull

   ; This call to cholesky_band operates on alpha and changes contents

   covariance = alpha
   errb = cholesky_band(alpha, mininf=min_influence) 
   if (errb[0] NE -1) then begin 
      if (arg_present(yfit)) then $
        yfit = bspline_valu(xdata, sset, x2=xdata, action=action, upper=upper, lower=lower)
      return, bspline_maskpoints(sset, errb, npoly)
   endif
 
   ; this changes beta to contain the solution

   errs = cholesky_solve(alpha, beta)   
   if (errs[0] NE -1) then begin
      if (arg_present(yfit)) then $
       yfit = bspline_valu(xdata, sset, x2=xdata, action=action, upper=upper, lower=lower)
      return, bspline_maskpoints(sset, errs, npoly)
   endif

   sc = size(sset.coeff)
   if (sc[0] EQ 2) then begin
      sset.icoeff[*,goodbk] = reform(alpha[0,lindgen(nfull)],npoly,nn)
      sset.coeff[*,goodbk] = reform(beta[lindgen(nfull)], npoly, nn)
   endif else begin
      sset.icoeff[goodbk] = alpha[0,lindgen(nfull)]
      sset.coeff[goodbk] = beta[lindgen(nfull)]
   endelse

   if (arg_present(yfit)) then $
    yfit = bspline_valu(xdata, sset, x2=xdata, action=action, upper=upper, lower=lower)

   return, 0L
end
