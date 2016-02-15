;+
; NAME:
;   mqladvance
;
; PURPOSE:
;   Generate MODIFIED associated Legendre polynomials one l at a time
;
; CALLING SEQUENCE:
;   mqladvance, x, m, Mql_1, Mql, lmax=lmax
;
; INPUTS:
;   x       - argument of Plm; in some cases x=cos(theta)
;   m       - m (scalar, double precision)
;   Mql_1   - Plm array, (Nx, lmax) of P_lm(x,l) for given m-1
;   
; KEYWORDS:
;   lmax    - set to lmax on first call; default determined from Mql_1
;
; OUTPUTS:
;   Mql     - Plm array, (Nx, lmax) of P_lm(x,l) for given m
;
; RESTRICTIONS:
;   Do not exceed m=lmax!
;
; EXAMPLES:
;   see healpix2alm.pro
;
; COMMENTS:
;   Based on prescription in Numerical Recipes. 
;
;   Must set either lmax (first time) or mql_1 (afterwards)
;
;   NOTE:  The Plms are the associated Legendre polynomials times 
;      sqrt((l-m)!/(l+m)!) for convenience in generating Ylms. 
;
;   First written years ago at Berkeley
;   Appears to be good to roughly machine roundoff error
;
;   This code is the same recursion as plmadvance, but runs
;     faster.  plmadvance is deprecated. 
;
; REVISION HISTORY:
;   2003-Feb-19  Written by Douglas Finkbeiner, Princeton
;   2003-Nov-13  Comments fixed up - DPF
;   2004-Aug-16  trap floating underflow - DPF
;
;----------------------------------------------------------------------
PRO mqladvance, x, m, Mql_1, Mql, lmax=lmax

; -------- Error checking
  if NOT keyword_set(lmax) then lmax = (size(Mql, /dimens))[1]
  if m gt lmax then begin 
     print, 'Your m value is too big!!!'
     return
  endif 

  if m eq 0 then begin 
     Mql = dblarr(n_elements(x), lmax+1)
     Mql[*, 0] = 1.d
     Mql[*, 1] = x
  endif else begin 

; -------- Seed recursion

; --------  m = l case, using P(l-1,m-1)
; NOTE: follow Numerical Recipes convention for minus sign;
;   IDL legendre function also contains this minus sign. 
;   Arfken does NOT have the minus. 

     lind = m
     l = double(m)
     Mql = Mql_1*0
     Mql[*, lind] = sqrt(1./((2*l)*(2*l-1))) * $
       (-1)*(2*l-1)*Mql_1[*, lind-1]*sqrt(1-x*x)
     cm = check_math()
     if (cm ne 0) and (cm ne 32) then stop
     if m eq lmax then return

  endelse 

; -------- Loop over l for given m, using P(l-1,m) and P(l-2,m)

  for lind=(m+1L) > 2L, lmax do begin 
     l = double(lind)
     temp = (2*l-1)*x*Mql[*, lind-1] - $
         ((l+m-1)*sqrt((l-m-1)/(l+m-1)))*Mql[*, lind-2] 
     Mql[*, lind] = sqrt((l-m)/(l+m))/(l-m) * temp
  endfor 

  cm = check_math()
  if (cm ne 0) and (cm ne 32) then stop

  return
END
