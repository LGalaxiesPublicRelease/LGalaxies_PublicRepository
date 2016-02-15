;+
; NAME:
;   hogg_iter_linfit
; PURPOSE:
;   Perform iterated linear least-squares fit with sigma-clipping.
; COMMENTS:
;   Solves the over-constrained equation yy=aa##xx, where the yy are
;     the "data" and the xx are the "parameters".
; CALLING SEQUENCE:
;   hogg_iter_linfit, aa,yy,ww,xx [,covar=covar]
; INPUTS:
;   aa      - input matrix, either transposition is fine
;   yy      - input vector of data in equation yy=aa##xx
;   ww      - weights for the components of yy; should be proportional to
;             (1/err^2); these are assumed to be exactly (1/err)^2 if keyword
;             "truesigma" is set
; OPTIONAL INPUTS:
;   nsigma  - number of sigma for clipping; defaults to 3.0
;   maxiter - maximum number of clipping iterations; defaults to 100
; KEYWORDS
;   median  - clip by deviation from median rather than from mean
;   verbose - output blow-by-blow
;   sacred  - don't update ww values to indicate clipped-out data points
;   truesigma - don't rescale weights and covariance matrix -- assume that the
;               supplied weights are exactly 1/sigma^2
; OUTPUTS:
;   xx      - output parameters in equation yy=aa##xx
;   covar   - covariance matrix for the fit parameters
;   ww      - ww values set to zero where data points clipped, and re-
;             scaled so that chisq^2 is unity per degree of freedom, unless
;             keyword "sacred" is set
; EXAMPLE:
;   Linear regression of ydata (with inverse variance ydata_ivar) on
;   xdata, with 5-sigma clipping:
;     aa= dblarr(2,n_elements(xdata))
;     aa[0,*]=1.D0
;     aa[1,*]=xdata
;     hogg_iter_linfit, aa, ydata, ydata_ivar, coeffs, covar=covar, nsigma=5
;   Resulting coefficients are in coeffs[0] and coeffs[1] 
; BUGS:
;   Covariance matrix is an approximation which assumes that the rms of the
;     fit is consistent with the errors in the yy, or something like that.
;     This could be changed if the input weights were always the true
;     1/sigma^2 values.
; REVISION HISTORY:
;   2000-Jun-03  written by Hogg, IAS
;   2000-Jun-12  output covariance matrix - Hogg
;   2000-Sep-06  added sacred keyword and updated ww matrix - Hogg
;   2001-Apr-10  allow rejected points to return to fit - Hogg
;-
;------------------------------------------------------------------------------
pro hogg_iter_linfit, aa,yy,ww,xx,covar=covar,nsigma=nsigma,maxiter=maxiter, $
                      median=median,verbose=verbose,sacred=sacred, $
                      truesigma=truesigma, condition=condition

; make everything double
   aa= double(aa)
   yy= double(yy)
   www= double(ww)

; check sizes and transposition
   if (size(aa))[0] EQ 1 then begin
     printf,-2, 'hogg_iter_linfit: converting aa to matrix form'
     aa= reform(aa,1,n_elements(aa))
   endif
   mm= (size(aa))[1]
   nn= (size(aa))[2]
   if nn EQ mm then begin
     printf,-2, 'hogg_iter_linfit: WARNING: aa matrix square'
     printf,-2, 'hogg_iter_linfit: WARNING: setting maxiter = 0'
     maxiter= 0
     printf,-2, 'hogg_iter_linfit: WARNING: and setting /truesigma'
     truesigma= 1
;     printf,-2, 'hogg_iter_linfit: ERROR: aa matrix square'
;     xx= sqrt(-1)
;     return
   endif
   if mm GT nn then begin
     printf,-2, 'hogg_iter_linfit: nn>mm, transposing aa matrix'
     aa= transpose(aa)
     mm= (size(aa))[1]
     nn= (size(aa))[2]
   endif
   if ((size(yy))[0] NE 1) OR ((size(yy))[1] NE nn) then begin
     printf,-2, 'hogg_iter_linfit: ERROR: yy vector incorrect size or shape'
     xx= sqrt(-1)
     return
   endif

; set iteration parameters
   if NOT keyword_set(nsigma) then nsigma= 3.0
   if NOT keyword_set(maxiter) then maxiter= 100
   isgood= byte(0*www+1)
   good= where(isgood,ngood)

; begin iterating
   iteration= 0
   repeat begin
     iteration= iteration+1

; do the matrix manipulations
     aat= transpose(aa)
     aatww= aat*((www*isgood)#(dblarr(mm)+1.0))
     aatwwaa= aatww##aa
     aatwwaainv= invert(aatwwaa)
     if mm GT 1 then aatwwyy= aatww##yy else aatwwyy= transpose(aatww)#yy
     if mm GT 1 then xx= aatwwaainv##aatwwyy else xx=aatwwaainv*aatwwyy

; make residuals
     residual= (yy-aa##xx)

; compute statistics
     if keyword_set(verbose) then $
       printf,-2, 'hogg_iter_linfit: iteration',iteration,': using', $
       ngood,' data points'
     if ngood GT 1 then begin
       if keyword_set(median) then begin
         medresid= median(residual[good])
         ms= mean((residual[good]-medresid)^2*www[good])
       endif else begin
         ms= mean((residual[good])^2*www[good])
       endelse
     endif else begin
       medresid= 0.0
       ms= 0.0
     endelse
     if keyword_set(verbose) then $
       printf,-2, 'hogg_iter_linfit: residual (weighted) sigma',sqrt(ms)

; sigma-clip and iterate
     oldngood= ngood
     if keyword_set(median) then begin
       isgood=((residual-medresid)^2*www LT nsigma^2*ms)
     endif else begin
       isgood=(residual^2*www LT nsigma^2*ms)
     endelse
     good= where(isgood,ngood)
     bad= where(isgood-1,nbad)
     if keyword_set(verbose) then $
       printf,-2, 'hogg_iter_linfit: will reject',nbad,' points'

; re-scale internal weights by rms (taking into account the number of
;   degrees of freedom) and iterate
     if NOT (keyword_set(truesigma) OR keyword_set(sacred)) then begin
       ndof= ngood-mm
       if ndof GT 0 then www= www*ndof/(ms*ngood)
     endif
   endrep until(((oldngood-ngood) EQ 0) OR ((nn-nbad) LE mm) $
                OR (iteration GE maxiter) OR (ngood EQ 0))

; make covariance matrix
   if keyword_set(truesigma) then begin
     covar= aatwwaainv
   endif else begin
     covar= aatwwaainv*ms*ngood/ndof
   endelse

; condition number of matrix
   if arg_present(condition) then begin
      if mm GT 1 then $
        condition = cond(aatwwaainv, /double, lnorm=2) $
      else condition = 1.0
   endif
   
; catch ugliness
   if (iteration GT 1) AND ((nn-nbad) LE mm) then begin
     printf,-2, 'hogg_iter_linfit: ERROR: sigma-clipped away too much data'
     xx= dblarr(mm)
     return
   endif

; update weights to internal values, if weights are not sacred
   if NOT keyword_set(sacred) then ww= www*isgood

   return
end
