;+
; NAME:
;	 nonneg_mult_update
; PURPOSE: (one line)
;	 Apply an SSL multiplicative update to iterate nonnegative quadratic problem 
; DESCRIPTION:
;  Using the method of Sha, Saul, & Lee (2002), "Multiplicative
;  updates for nonnegative quadratic programming in support vector
;  machines" (UPenn Tech Report MS-CIS-02-19), apply a multiplicative
;  update to an attempted solution to a nonnegative quadratic problem
;  (QP with box constraints):
;
;     F(v) = (1/2) v^T.A.v + b^T.v for v_i >= 0 for all i.
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
; CATEGORY:
;       Numerical
; CALLING SEQUENCE:
;	 new = nonneg_mult_update(old,avfunc,b)
; INPUTS:
;	 old - start vector
;  avfunc - function which returns A+.v or A-.v, depending
;  b - vector
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;  factor - return the factor used in this vector
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
;  Untested.
; PROCEDURE:
; MODIFICATION HISTORY:
;	  Written 2003-01-02 MRB (NYU) at suggestion of Sam Roweis
;-
pro nonneg_mult_update,old,new,avfunc,b,factor=factor

if(n_params() ne 4) then begin
    print,'Syntax - new=nonneg_mult_update(old, avfunc, b [, factor=])'
    return
endif 

avpos=call_function(avfunc,old,1.)
avneg=call_function(avfunc,old,-1.)

; if you are at zero, multiplicative updates don't change
; your value, so just ignore
nnindx=where(old gt 0.D,nncount)
if(nncount gt 0) then $
    new[nnindx]=old[nnindx]* $
      (-b[nnindx]+sqrt(b[nnindx]^2+4.*avpos[nnindx]*avneg[nnindx]))/ $
      (2.*avpos[nnindx])

return

end
