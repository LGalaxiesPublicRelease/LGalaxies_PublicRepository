;+
; NAME:
;     adjust_error
; PURPOSE:
;     given data points, uncertainties, and mean value, determines
;     what extra uncertainty yields chi^2/dof=1?
; CALLING SEQUENCE
;      adjustment = adjust_error( x, xerr, mean=mean, dof=dof)
; INPUTS:
;     x - [N] array of measurements
;     xerr - [N] array of estimated measurement uncertainties
; OPTIONAL INPUTS:
;     mean - mean value if known [default is computed from x, xerr]
;     dof - degrees of freedom, defaults to N; value (even the
;           default) is reduced by 1 if mean is not input
;     tol - tolerance in adjustment in chi^2 units [default is
;           0.01*sqrt(dof)]
;     maxiter - maximum number of binary search iterations
; OUTPUT:
;     adjustment - minimum uncertainty to add to each xerr (in
;                  quadrature) to yield chi^2 <= dof 
; OPTIONAL OUTPUTS:
;     niter - number of iterations executed
;     chisq - chi^2 corresponding to adjustment
; METHOD:
;     dunno
; REVISION HISTORY:
;     Blanton & Hogg 2003-07-15 written, tested a bit
;-
function calc_chisq, x, xerr, mean, adjustment
return,total((x-mean)^2/(xerr^2+adjustment^2),/double)
end
;
function adjust_error, x, xerr, mean=mean, dof=dof, tol=tol, niter=niter, $
                       maxiter=maxiter, chisq=chisq

if(not keyword_set(maxiter)) then  maxiter=10000000L
nx=n_elements(x)
if(n_elements(dof) eq 0) then dof=nx
ivar=1.D/xerr^2
if(n_elements(mean) eq 0) then begin
    mean=total(x*ivar,/double)/total(ivar,/double)
    dof=dof-1L
endif
if(not keyword_set(tol)) then $
  tol=0.01/sqrt(double(dof))

; initial to bracket the root
niter=0
adjustment0=0.D
chisq=calc_chisq(x,xerr,mean,adjustment0)
if(chisq lt dof) then return,adjustment0
adjustment1=median(xerr)
chisq1=calc_chisq(x,xerr,mean,adjustment1)
while(chisq1 ge dof) do begin
    adjustment0=adjustment1
    adjustment1=2.D*adjustment1
    chisq1=calc_chisq(x,xerr,mean,adjustment1)
endwhile

while(abs(chisq-double(dof)) gt tol AND niter lt maxiter) do begin
    adjustmentmid=0.5D*(adjustment0+adjustment1)
    chisq=calc_chisq(x,xerr,mean,adjustmentmid)
    if(chisq lt dof) then $
      adjustment1=adjustmentmid $
    else $
      adjustment0=adjustmentmid 
    niter=niter+1L
endwhile

return, adjustmentmid
  
end
