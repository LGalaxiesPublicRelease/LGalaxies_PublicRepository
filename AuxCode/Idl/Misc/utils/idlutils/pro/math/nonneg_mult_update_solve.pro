;+
; NAME:
;	 nonneg_mult_update_solve
; PURPOSE: (one line)
;	 Use nonneg_mult_update to iterate to convergence
; DESCRIPTION:
;  From some starting point, iterates to convergence a
;  box-constrained QP problem: 
;
;     F(v) = (1/2) v^T.A.v + b^T.v for v_i >= 0 for all i.
;
;  Uses the method of Sha, Saul, & Lee (2002), "Multiplicative
;  updates for nonnegative quadratic programming in support vector
;  machines" (UPenn Tech Report MS-CIS-02-19).
;
;  It requires the user to supply a function avfunc(vec,sign) which returns
;  A+.v if sign>0. and A-.v if sign<0, where:
;  
;     A+_ij = A_ij for A_ij>0.
;             0.   otherwise
;
;     A-_ij = |A_ij|  for A_ij<0.
;             0.      otherwise
; 
;  Alternatively, if /matrix is set, nonneg_mult_update_solve will
;  interpret the input avfunc as a matrix and construct the
;  appropriate functions. 
; CATEGORY:
;       Numerical
; CALLING SEQUENCE:
;	 solution = nonneg_mult_update_solve(start,avfunc,b)
; INPUTS:
;	 start - start vector
;  avfunc - function which returns A+.v or A-.v, depending
;  b - vector
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
;  /matrix - indicates that avfunc is actually just the matrix A,
;            and the default A+/-.v function should be used
; OUTPUTS:
;	 Return value is the shifted array.
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
;  Tested only in simple cases.
; PROCEDURE:
; MODIFICATION HISTORY:
;	  Written 2003-01-02 MRB (NYU) at suggestion of Sam Roweis
;-
function nnmus_default_avfunc,vec,sign

common nnmus_com, nnmus_matrix_pos, nnmus_matrix_neg

if(sign gt 0.) then return, nnmus_matrix_pos#vec
return, nnmus_matrix_neg#vec

end
;
function nnmus_value,vec,avfunc,b,offset=offset

common nnmus_com, nnmus_matrix_pos, nnmus_matrix_neg

avpos=call_function(avfunc,vec,1.)
avneg=call_function(avfunc,vec,-1.)
av=avpos-avneg
vav=(transpose(vec)#av)[0]

val=(0.5*vav+transpose(b)#vec)[0]
if(keyword_set(offset)) then val=val+offset
return,val

end
;
function nonneg_mult_update_solve,start,avfunc,b,matrix=matrix,tol=tol, $
                                  verbose=verbose,offset=offset,value=value, $
                                  maxiter=maxiter,niter=niter,skip=skip, $
                                  chi2tol=chi2tol

common nnmus_com

if(n_params() ne 3) then begin
    print,'Syntax - result=nonneg_mult_update_solve(start, avfunc, b [, /matrix, tol=, $'
    print,'             /verbose, offset=, value=, maxiter=, niter=])'
    return,-1
endif

if(NOT keyword_set(skip)) then skip=10
if(NOT keyword_set(tol)) then tol=1.D-7
b=reform(b,n_elements(b))
start=reform(start,n_elements(start))

; set avfunc to use
use_avfunc=avfunc
if(keyword_set(matrix)) then begin
    nnmus_matrix_pos=avfunc > 0.D
    nnmus_matrix_neg=abs(avfunc < 0.D)
    use_avfunc='nnmus_default_avfunc'
endif

sol=start
oldval=nnmus_value(sol,use_avfunc,b,offset=offset)
if(keyword_set(verbose)) then splog,'oldval='+string(oldval)
diff=tol*2.
niter=0
next=dblarr(n_elements(sol))
while(abs(diff) gt tol) do begin
    nonneg_mult_update,sol,next,use_avfunc,b,factor=factor
    sol=next
    if((niter mod skip) eq 0) then begin
        newval=nnmus_value(sol,use_avfunc,b,offset=offset)
        diff=(newval-oldval)
        oldval=newval
    endif
    ;if(keyword_set(verbose)) then begin
        ;splog,'newval='+string(newval)
;        splog,'sol='+string(sol)
        ;splog,['min(factor)=','max(factor)']+string(minmax(factor))
    ;endif
    niter=niter+1L
    if(keyword_set(maxiter)) then $
      if(niter ge maxiter) then $
      diff=0.D
    if(keyword_set(chi2tol)) then $
       if(oldval lt chi2tol) then $
       diff=0.D
endwhile

value=oldval

return,sol

end
