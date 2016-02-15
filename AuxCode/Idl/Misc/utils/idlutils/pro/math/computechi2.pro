;+
; NAME:
;   computechi2
;
; PURPOSE:
;   Solve the linear set of equations Ax=b using SVD
;
; CALLING SEQUENCE:
;   chi2 = computechi2( bvec, sqivar, amatrix, $
;    [ acoeff=, dof=, yfit=, covar=, var= ] )
;
; INPUTS:
;    bvec      - b vector in Ax=b [N]
;    sqivar    - Errors in b as 1/sigma [N]
;    amatrix   - A matrix in Ax=b [N,M]
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   chi2       - Chi^2 of the fit
;
; OPTIONAL OUTPUTS:
;   acoeff     - Fit parameters x in Ax=b [M]
;   dof        - Degrees of freedom in the fit, equal to the number of
;                equations where SQIVAR is not zero minus the number of
;                fit parameters (M)
;   yfit       - Evaluation of the best-fit at each data point [N]
;   covar      - Covariance matrix [M,M]
;   var        - Variances [M], which is equivalent to the diagonal
;                entries of COVAR
;
; COMMENTS:
;
; EXAMPLES:
;   The following example creates a data vector with 20 measurements
;   and their errors.  A linear fit is performed using the LINFIT()
;   function and this function, which are equivalent:
;     n = 20
;     x = dindgen(n)
;     y = 10.d0 * smooth(randomu(-1234,n),10)
;     sqivar = 0.5 + randomu(-4321,n)
;     acoeff = linfit(x,y,measure_errors=1./sqivar,covar=covar)
;     print, acoeff, covar
;     templates = [[dblarr(n)+1],[dindgen(n)]]
;     chi2 = computechi2(y,sqivar,templates,acoeff=acoeff,covar=covar)
;     print, acoeff, covar
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   07-Aug-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
function computechi2, objflux, sqivar, starflux, $
 acoeff=acoeff, dof=dof, yfit=yfit, covar=covar, var=var

   ndim = size(starflux, /n_dimen)
   if (ndim EQ 1) then nstar = 1 $
    else nstar = (size(starflux, /dimens))[1]

   bvec = double(objflux) * double(sqivar)
   mmatrix = double(starflux) * (double(sqivar) # replicate(1,nstar))

;   ---------  the line above is about twice as fast --------------
;   mmatrix = starflux
;   for i=0L, nstar-1 do $
;     mmatrix[*,i] = mmatrix[*,i] * sqivar

   mmatrixt = transpose( mmatrix )
   mm = mmatrixt # mmatrix

   ; Use SVD to invert the matrix
;   mmi = invert(mm, /double)
   if (nstar EQ 1) then begin
      ; The last term below is to protect against divide-by-zero
      ; in the degenerate case.
      mmi = 1.0 / (mm + (mm EQ 0))
   endif else begin
      svdc, mm, ww, uu, vv, /double
      mmi = 0 * vv
      ; The last term below is to protect against divide-by-zero
      ; in the degenerate case.
      for i=0L, nstar-1 do mmi[i,*] = vv[i,*] / (ww[i] + (ww[i] EQ 0))
      mmi = mmi ## transpose(uu)
   endelse

   acoeff = mmi # (mmatrixt # bvec)
   chi2 = total( (mmatrix # acoeff - bvec)^2, /double )

   if (arg_present(yfit)) then $
    yfit = acoeff ## double(starflux)
;    yfit = transpose(acoeff # mmatrixt) / sqivar
   if (arg_present(dof)) then $
    dof = total(sqivar NE 0) - nstar

;   if (arg_present(covar)) then begin
;      covar = dblarr(nstar,nstar)
;      for j=0L, nstar-1L do begin
;         for k=0L, nstar-1L do begin
;            covar[j,k] = total(vv[j,*]^2 * vv[k,*]^2 / (ww[k]), /double)
;         endfor
;      endfor
;   endif

   if (nstar EQ 1) then begin
      ivar = sqivar^2
      wtot = total(ivar,/double)
      if (wtot LE 0) then begin
         var = [0]
      endif else begin
         xmean = total(objflux * ivar, /double) / wtot
         var = [total((objflux - xmean)^2 * ivar, /double)]
      endelse
      covar = var
   endif else begin
      if (arg_present(covar) OR arg_present(var)) then begin
         igood = where(ww GT 0, ngood)
         wwt = dblarr(nstar)
         if (ngood GT 0) then wwt[igood] = 1.d0 / ww[igood]
      endif

      if (arg_present(covar)) then begin
         covar = dblarr(nstar,nstar)
         for i=0L, nstar-1L do begin
            for j=0L, i do begin
               covar[i,j] = total(wwt * vv[*,i] * vv[*,j], /double)
               covar[j,i] = covar[i,j]
            endfor
         endfor
      endif

      if (arg_present(var)) then begin
         var = dblarr(nstar)
         if (arg_present(covar)) then begin
            i = lindgen(nstar)
            var[*] = covar[i,i]
         endif else begin
            for j=0L, nstar-1L do $
             var[j] = total((vv[j,*])^2 * wwt, /double)
         endelse
      endif
   endelse

   return, chi2
end
;------------------------------------------------------------------------------
