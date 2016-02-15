;+
; NAME:
;	  gauss_overlap
; PURPOSE: (one line)
;   calculate log integral, covariance, and mean of product of two gaussians
; DESCRIPTION:
;   calculate log integral, covariance, and mean of product of two gaussians
; CATEGORY:
;       Numerical
; CALLING SEQUENCE:
;	  gauss_overlap,mean1,covar1,mean2,covar2,mean=mean,covar=covar, $
;	    logint=logint, invcovar=invcovar, /input_invcovar
; INPUTS:
;	  mean1 - [N] mean of first gaussian
;	  covar1 - [N, N] covariance of first gaussian
;	  mean2 - [N] mean of second gaussian
;	  covar2 - [N, N] covariance of second gaussian
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
;   /input_invcovar - covar1, covar2 are actually inverse covars
; OUTPUTS:
;   mean - mean of global gaussian
;   covar - covariance of global gaussian
;   invcovar - inverse covariance of global gaussian
;   logint - natural log of integral of global gaussian
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; MODIFICATION HISTORY:
;	  Blanton and Roweis 2003-02-18j 
;-
pro gauss_overlap,mean1,covar1,mean2,covar2,mean=mean,covar=covar, $
                  invcovar=invcovar, logint=logint, $
                  input_invcovar=input_invcovar

ndim=n_elements(mean1)

if(keyword_set(invcovar)) then begin
    invcovar1=double(covar1)
    invcovar2=double(covar2)
endif else begin
    invcovar1=invert(double(covar1),/double)
    invcovar2=invert(double(covar2),/double)
endelse
invcovar=invcovar1+invcovar2
covar=invert(invcovar,/double)
mean=covar#(invcovar1#double(mean1)+invcovar2#double(mean2))

if(ndim eq 1) then begin
    minus_log_z1=0.5*alog(invcovar1)
    minus_log_z2=0.5*alog(invcovar2)
    minus_log_z=0.5*alog(invcovar)
endif else begin
    minus_log_z1=0.5*alog(determ(invcovar1))
    minus_log_z2=0.5*alog(determ(invcovar2))
    minus_log_z=0.5*alog(determ(invcovar))
endelse
muCmu1=double(mean1)#invcovar1#double(mean1)
muCmu2=double(mean2)#invcovar2#double(mean2)
muCmu=mean#invcovar#mean

logint=-0.5*double(ndim)*alog(2.*!DPI)+ $
  minus_log_z1+ $
  minus_log_z2- $
  minus_log_z- $
  muCmu1- $
  muCmu2+ $
  muCmu

end
