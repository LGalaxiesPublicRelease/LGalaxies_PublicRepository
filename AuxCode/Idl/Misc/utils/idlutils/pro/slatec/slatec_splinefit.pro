;+
; NAME:
;   slatec_splinefit
;
; PURPOSE:
;   Calculate a B-spline in the least squares sense with rejection
;
; CALLING SEQUENCE:
;   
;   fullbkpt = slatec_splinefit(x, y, coeff, invvar=, x2=, npoly=, $
;    bkpt=, fullbkpt=, maxiter=, upper=, lower=, $
;    rejper=, /eachgroup, /secondkludge, $
;    mask=, _EXTRA=KeywordsForEfc)
;
; INPUTS:
;   x          - Data x values
;   y          - Data y values
;   bkpt       - Breakpoint vector returned by SLATEC_EFC()
;
; OPTIONAL KEYWORDS:
;   invvar     - Inverse variance of y; if not set, then set to be
;                consistent with the standard deviation.  This only matters
;                if rejection is being done.
;   x2         - 2nd dependent variable for 2-D fitting 
;   npoly      - Polynomial order to fit over 2nd variable (X2); default to 2.
;   maxiter    - maximum number of iterations (default 5)
;   upper      - Sigma rejection threshold for positive deviations; default to 5
;   lower      - Sigma rejection threshold for negative deviations; default to 5
;   rejper     - Alternative rejection algorithm, rejecting at most the
;                fraction REJPER of the points per iteration (but rejecting
;                at least one point if there are any bad ones).
;   eachgroup  - Alternative rejection algorithm ???
;   secondkludge - ???
;
; KEYWORDS FOR SLATEC_EFC:
;   nord
;   bkspace
;   nbkpts
;   everyn
;
; OUTPUTS:
;   fullbkpt   - The fullbkpt vector required by evaluations with
;                SLATEC_BVALU().  Return -1 if the spline fit fails
;                due to all points (or all but 1) being rejected.
;
; OPTIONAL OUTPUTS:
;   coeff      - B-spline coefficients calculated by SLATEC_EFC().
;   bkpt       - Breakpoints without padding
;   mask       - Mask array, set to 0 for good points and 1 for rejected points
;
; COMMENTS:
;   If both bkspace and nbkpts are passed, bkspace is used
;
; EXAMPLES:
;
;   x = findgen(100)
;   y = randomu(100,100)
;   fullbkpt = slatec_splinefit(x, y, coeff, invvar=invvar, nbkpts=10)
;
;   xfit = findgen(10)*10.0
;   yfit = slatec_bvalu(xfit, fullbkpt, coeff)
;
; PROCEDURES CALLED:
;   slatec_bvalu()
;   slatec_efc()
;
; REVISION HISTORY:
;   15-Oct-1999  Written by Scott Burles, Chicago
;-
;------------------------------------------------------------------------------
function slatec_splinefit, x, y, coeff, invvar=invvar, x2=x2, npoly=npoly, $
 bkpt=bkpt, fullbkpt=fullbkpt, maxiter=maxiter, upper=upper, lower=lower, $
 rejper=rejper, eachgroup=eachgroup, secondkludge=secondkludge, $
 mask=mask, _EXTRA=KeywordsForEfc

   if (n_params() LT 3) then begin
      print, 'Syntax -  fullbkpt = slatec_splinefit(x, y, coeff, invvar=, x2=, npoly=, $'
      print, '           bkpt=, fullbkpt=, maxiter=, upper=, lower=, $'
      print, '           rejper=, /eachgroup, /secondkludge, $'
      print, '           mask=, _EXTRA=KeywordsForEfc)'
   endif

   if (n_elements(x) NE N_elements(y)) then $
    message, 'Dimensions of X and Y do not agree'
   if (n_elements(maxiter) EQ 0) then maxiter = 5
   if (NOT keyword_set(upper)) then upper = 5.0
   if (NOT keyword_set(lower)) then lower = 5.0
   if (NOT keyword_set(invvar)) then begin
      ; Set INVVAR to be consistent with the standard deviation
      invvar = y - y + 1.0 / (stddev(y,/double))^2
   endif else begin
      if (N_elements(y) NE N_elements(invvar)) then $
       message, 'Dimensions of Y and INVVAR do not agree'
   endelse

   invsig = sqrt(invvar > 0)

   nx = n_elements(x)
   mask = bytarr(nx) + 1
   bad = where(invvar LE 0.0)
   good = where(invvar GT 0.0)
   if (bad[0] NE -1) then mask[bad] = 0

   for iiter=0, maxiter do begin
      oldmask = mask
      these = where(mask, nthese)

      if (nthese LT 2) then begin
         print, 'Lost all points'
         return, -1
      endif

      if (NOT keyword_set(x2)) then begin
         fullbkpt = slatec_efc(x[these], y[these], fullbkpt=fullbkpt, $
              coeff, bkpt=bkpt, invsig=invsig[these], $
                 _EXTRA=KeywordsForEfc)
         yfit = slatec_bvalu(x, fullbkpt, coeff)
      endif else begin
         fullbkpt = slatec_efc(x[these], y[these], fullbkpt=fullbkpt, $
              coeff, bkpt=bkpt, invsig=invsig[these], x2=x2[these], $
              npoly=npoly, _EXTRA=KeywordsForEfc)
         yfit = bvalu2d(x, x2, fullbkpt, coeff)
      endelse

      if (keyword_set(rejper)) then begin
         diff = (y[these] - yfit[these])*invsig[these]
         bad = where(diff LE -lower OR diff GE upper, nbad)
         if (nbad EQ 0) then iiter = maxiter $
         else begin
            negs = where(diff LT 0)
            tempdiff = abs(diff)/upper
            if (negs[0] NE -1) then tempdiff[negs] = abs(diff[negs])/lower

            worst = reverse(sort(tempdiff))
            rejectspot = long(nbad * rejper)
            mask[these[worst[0:rejectspot]]] = 0

            good = where(diff GT -lower AND diff LT upper AND $
                        invsig[these] GT 0.0)
            if (good[0] NE -1) then mask[these[good]] = 1
         endelse
      endif else if (keyword_set(eachgroup)) then begin
         diff = (y[these] - yfit[these])*invsig[these]
         bad = where(diff LE -lower OR diff GE upper, nbad)
         if (nbad EQ 0) then iiter = maxiter $
         else begin
            negs = where(diff[bad] LT 0)
            tempdiff = abs(diff[bad])/upper
            if (negs[0] NE -1) then tempdiff[negs] = abs(diff[bad[negs]])/lower
            if (nbad EQ 1) then mask[these[bad]] = 0 $
            else begin
               groups = where(bad[1:nbad-1] - bad[0:nbad-2] NE 1, ngroups)
               if (ngroups EQ 0) then groups = [-1,nbad-1] $
               else groups = [-1,groups,nbad-1]
               for i=0L, ngroups do begin
                  groupmax = max(tempdiff[groups[i]+1:groups[i+1]], place)
                  mask[these[bad[place+groups[i]+1]]] = 0
               endfor
            endelse 

            good = where(diff GT -lower AND diff LT upper AND $
                        invsig[these] GT 0.0)
            if (good[0] NE -1) then mask[these[good]] = 1
         endelse
      endif else begin
         diff = (y - yfit)*invsig
         bad = where(diff LT -lower OR diff GT upper OR invsig LE 0.0, nbad)
         if (nbad EQ 0) then iiter = maxiter $
         else begin
            mask = bytarr(nx) + 1
            mask[bad]  = 0
            if total(abs(mask - oldmask)) EQ 0 then iiter=maxiter
         endelse

         ;----------
         ; End iterations if no points are rejected, or the same points are
         ; rejected as the last iteration.

         if (nbad EQ 0) then iiter = maxiter $
         else begin
            mask = bytarr(nx) + 1
            mask[bad]  = 0
            if total(abs(mask - oldmask)) EQ 0 then iiter=maxiter
         endelse
      endelse

   endfor

   mask = oldmask ; Return the mask that was actually used in the last
                  ; call to SLATEC_EFC()

   if keyword_set(secondkludge) then begin
      nfullbkpt = n_elements(fullbkpt) 
      ncoeff = n_elements(coeff) 
      nordtemp = nfullbkpt - ncoeff
      nbkptfix = nfullbkpt - 2*nordtemp + 1

      xbkpt = fullbkpt[nordtemp:nordtemp+nbkptfix-1]
      deriv2 = sqrt(abs(slatec_bvalu(x, fullbkpt, coeff)))
      deriv2total = deriv2
      for i=1L, nx-1 do $
       deriv2total[i] = deriv2total[i] + deriv2total[i-1]
      total2 = total(deriv2)
      deriv2x = (findgen(nbkptfix)+1.0)*total2/float(nbkptfix)
      invspl = spl_init(deriv2total,x)
      newbkpt = spl_interp(deriv2total,x,invspl,deriv2x)

      fullbkpt[nordtemp:nordtemp+nbkptfix-1] = newbkpt

      these = where(mask)
      fullbkpt = slatec_efc(x[these], y[these], fullbkpt=fullbkpt, $
       coeff, invsig=invsig[these], _EXTRA=KeywordsForEfc)
   endif

   return, fullbkpt
end
