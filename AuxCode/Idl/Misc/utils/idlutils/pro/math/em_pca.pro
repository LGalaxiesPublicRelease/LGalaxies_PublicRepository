;+
; NAME:
;    em_pca
; PURPOSE: (one line)
;    Perform E-M PCA to get first k principal components
; DESCRIPTION:
;    Uses Sam Roweis' Expectation Maximization version of PCA to
;    efficiently find the first k principal components for a
;    distribution in p (>k) dimensions. The procedure guesses an
;    initial set of eigenvectors (stored the in [p,k] dimensional
;    matrix "eigenvec") and applies the following iteration to the
;    [p,N] dimensional "data":
;
;         hidden= ( eigenvec^T . eigenvec )^{-1} . eigenvec^T . data
;         eigenvec= data . hidden^T . (hidden . hidden^T)^{-1}
;
;    From:
;    Neural Information Processing Systems 10 (NIPS'97) pp.626-632
;    available at:
;    http://www.cs.toronto.edu/~roweis/papers/empca.ps.gz
; CATEGORY:
;    Mathematical
; CALLING SEQUENCE:
;    em_pca, data, k, eigenvec, hidden [, tol=, maxiter=, niter=, /verbose]
; INPUTS:
;    data - [p,N] data to be PCAed
;    k - number of eigenvectors desired (<p)
; OPTIONAL INPUT PARAMETERS:
;    tol - tolerance of convergence (default 0.)
;    maxiter - maximum number of iterations (default 20)
; KEYWORD PARAMETERS:
;    /verbose - verbose output
;    /nofix - don't do the final real PCA
; OUTPUTS:
;    eigenvec - [p,k] matrix of k leading eigenvectors
;    hidden - [k, N] matrix of "hidden" variables (the lower dimensional
;             representation of the data)
; OPTIONAL OUTPUTS:
;    niter - number of iterations used
; COMMON BLOCKS:
; SIDE EFFECTS:
; BUGS:
;    Somewhat untested.
; RESTRICTIONS:
;    Does not implement Sam's Sensible-PCA procedure
; PROCEDURE:
; MODIFICATION HISTORY:
;    2003-01-26 - Written by Michael Blanton (NYU)
;-
pro em_pca, data, k, eigenvec, hidden, tol=tol, maxiter=maxiter, niter=niter, $
            verbose=verbose, nofix=nofix, noortho=noortho

; set defaults
if(n_elements(tol) eq 0) then tol=0.
if(n_elements(maxiter) eq 0) then maxiter=20

; check args
if (n_params() lt 1) then begin
    print, 'Syntax - em_pca, data, k, eigenvec, hidden [, tol=, maxiter=, niter=, $' 
    print, '                 /verbose]'
    return
endif

; check dimensions
ndim=(size(data))[0]
if(ndim eq 1) then $
  p=1 $
else $
  p=(size(data,/dimensions))[0]
n=n_elements(data)/p
if(k gt p) then begin
    splog, 'k must be less than or equal to p!'
    splog, 'k= '+string(k)
    splog, 'p= '+string(p)
    return
endif

; initial conditions
eigenvec=dblarr(p,k)
eigenvec[0:k-1,0:k-1]=identity(k)

niter=0L
diff=tol*2.+1.
while(niter lt maxiter and diff gt tol) do begin
    if(keyword_set(verbose)) then splog,'niter= '+string(niter)
    hidden=(invert(transpose(eigenvec)#eigenvec,/double)# $
            transpose(eigenvec))#data
    oldeigenvec=eigenvec
    eigenvec=data#transpose(hidden)#invert(hidden#transpose(hidden),/double)
    if(tol gt 0.) then begin
        diff=0.
        for i=0, k-1 do begin
            diff=diff+abs(1.-total(oldeigenvec[*,i]*eigenvec[*,i],/double)/ $
                          sqrt(total(oldeigenvec[*,i]^2,/double)* $
                               total(eigenvec[*,i]^2,/double)))
        endfor
        if(keyword_set(verbose)) then splog,'diff= '+string(diff)
    endif
    niter=niter+1L
endwhile

if(NOT keyword_set(noortho)) then begin
;   Orthonormalize
    for b = 0l, k-1l do begin
;       orthogonalize
        for bp = 0l, b-1l do begin
            dot=total(eigenvec[*,b]*eigenvec[*,bp],/double)
            eigenvec[*,b]=eigenvec[*,b]-dot*eigenvec[*,bp]
        endfor 
        
;       normalize
        dot=total(eigenvec[*,b]^2,/double)
        dot=1./sqrt(dot)
        eigenvec[*,b]=eigenvec[*,b]*dot
    endfor 

    if(NOT keyword_set(nofix)) then begin
;       project variables onto new coordinates
        hidden=transpose(eigenvec)#data
        pca, transpose(hidden), eval_hidden, evec_hidden, /silent, /covariance
        eigenvec=eigenvec[*,*]#transpose(evec_hidden)
    endif

;   project variables onto new coordinates
    hidden=transpose(eigenvec)#data
endif

return

end
