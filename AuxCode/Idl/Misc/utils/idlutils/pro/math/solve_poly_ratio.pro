;+
; NAME:
;   solve_poly_ratio
;
; PURPOSE:
;   Compute the polynomial function that multiplies one spectrum
;   to agree with another
;
; CALLING SEQUENCE:
;   solve_poly_ratio( xvector, aflux, bflux, aivar, [ bivar, npoly=, nback=, $
;    inparams=, yfit=, ymult=, yadd=, acoeff=, totchi2=, status=, perror ] )
;
; INPUTS:
;   xvector    - X axis from which to construct polynomial terms; this vector
;                should be well-conditioned (for example, in the range [0,1])
;                such that raising it to powers does not result in numeric
;                overflows [NPIX]
;   aflux      - Flux vector to be rescaled [NPIX]
;   bflux      - Reference flux vector [NPIX]
;   aivar      - Inverse variance for AFLUX [NPIX]
;
; OPTIONAL INPUTS:
;   bivar      - Inverse variance for BFLUX [NPIX]
;   npoly      - Number of polynomial terms for multiplying AFLUX; this must
;                be a positive integer; default to 1
;   nback      - Number of polynomial terms to add to AFLUX
;   inparams   - Starting guess for polynomial + additive terms if the
;                second method is used (specifying AIVAR,BIVAR) [NPOLY+NBACK];
;                if not set, then a fit is first performed without
;                allowing the errors AIVAR to be rescaled with the flux
;   status     - Return value from MPFIT if using the 2nd fit method
;   perror     - Return value from MPFIT if using the 2nd fit method
;
; OUTPUTS:
;   yfit       - Rescaled AFLUX with additive, background terms:
;                AFLUX * SUM{i=0,NPOLY-1} XVECTOR^i + SUM{j=0,NBACK-1} XVECTOR^j
;
; OPTIONAL OUTPUTS:
;   ymult      - Final multiplicative polynomial [NPIX]
;   yadd       - Final addditive polynomial [NPIX]
;   acoeff     - All coefficicents, starting with the multiplicative
;                [NMULT+NBACK]
;   totchi2    - Total chi^2
;
; COMMENTS:
;   Fit for the polynomial vector multiplying one vector to agree
;   with another such that
;     BFLUX = polynomial * AFLUX + background
;
;   The polynomials are defined such that:
;      YFIT = YMULT * AFLUX + YADD
;
;   There are two modes of operation for this routine.
;   The simpler algorithm is if AIVAR is specified but not BIVAR.
;   One is minimzing the sum of the squares of the chi values:
;     chi = (BFLUX - polynomial * AFLUX - background) * sqrt(AIVAR)
;   That is a linear problem that we solve with SVD.
;
;   The second algorithm has an error associates with each input vector.
;   The problem is now nonlinear, because BIVAR is rescaled appropriately
;   with BFLUX.  The errors are the quadrature sum from AIVAR and the
;   rescaled BIVAR.
;   The errors are rescaled at every function evaluation, but we
;   only allow the errors to get smaller by up to a factor of 1e4,
;   and we only allow them to get larger slowly (as the square root).
;   This should very strongly constrain the flux-corrrection vectors
;   from going too small (or negative), or too large.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   mpfit()
;
; INTERNAL SUPPORT ROUTINES:
;   solve_poly_fn()
;   solve_poly_chi_fn()
;   solve_poly_vectors()
;   solve_poly_ratio2()
;   solve_poly_ratio1()
;
; REVISION HISTORY:
;   05-Feb-2004  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
forward_function mpfit, solve_poly_chi_fn

;------------------------------------------------------------------------------
function solve_poly_fn, acoeff, ymult=ymult, yadd=yadd

   common com_fcorr_chi, npoly, nback, xvector, aflux, bflux, $
    aivar, bivar, aarr, barr

   ymult = acoeff[0] * aarr[*,0]
   for i=1, npoly-1 do ymult = ymult + acoeff[i] * aarr[*,i]
   if (nback EQ 0) then yadd = 0 $
    else yadd = acoeff[npoly:npoly+nback-1] ## barr
   yfit = ymult * aflux + yadd

   return, yfit
end
;------------------------------------------------------------------------------
; Return a vector of chi values
function solve_poly_chi_fn, acoeff

   common com_fcorr_chi, npoly, nback, xvector, aflux, bflux, $
    aivar, bivar, aarr, barr

   yfit = solve_poly_fn(acoeff, ymult=ymult)

   ; The errors are rescaled at every function evaluation, but we
   ; only allow the errors to get smaller by up to a factor of 1e4,
   ; and we only allow them to get larger slowly (as the square root).
   ;  This should very strongly constrain the flux-corrrection vectors
   ; from going too small (or negative), or too large.
   qgood = aivar GT 0 AND bivar GT 0
   vmult = (ymult > 1e-4) * (ymult LE 1) + sqrt(ymult) * (ymult GT 1)
   totivar = qgood / ( 1./(aivar + (1-qgood)) + vmult^2/(bivar + (1-qgood)) )

   chivec = (bflux - yfit) * sqrt(totivar)

   return, chivec
end

;------------------------------------------------------------------------------
; Construct the polynomial or additive vectors
function solve_poly_vectors, xvector, npoly

   if (npoly EQ 0) then return, 0

   npix = n_elements(xvector)
   aarr = dblarr(npix, npoly)
   for ipoly=0, npoly-1 do aarr[*,ipoly] = xvector^ipoly

   return, aarr
end
;------------------------------------------------------------------------------
pro solve_poly_ratio2, xvector1, allflux1, allflux2, allivar1, allivar2, $
 npoly=npoly1, nback=nback1, yfit=yfit, ymult=ymult, yadd=yadd, $
 acoeff=acoeff, totchi2=totchi2, $
 inparams=inparams, status=status, perror=perror

   common com_fcorr_chi, npoly, nback, xvector, aflux, bflux, $
    aivar, bivar, aarr, barr

   ; Pass arrays in the common block
   if (keyword_set(npoly1)) then npoly = npoly1[0] $
    else npoly = 1
   if (n_elements(nback1) GT 0) then nback = nback1[0] $
    else nback = 0
   xvector = xvector1
   aflux = allflux1
   bflux = allflux2
   aivar = allivar1
   bivar = allivar2

   aarr = solve_poly_vectors(xvector, npoly)
   barr = solve_poly_vectors(xvector, nback)

   ; Call MPFIT to iterate on the solution for the template
   parinfo1 = {value: 0.D, fixed: 0, limited: [0b,0b], limits: [0.d0,0.d0], $
    tied: ''}
   parinfo = replicate(parinfo1, npoly+nback)
   if (keyword_set(inparams)) then begin
      parinfo[0:n_elements(inparams)-1].value = inparams
   endif else begin
      parinfo[0].value = 1.D
      if (npoly GT 1) then parinfo[1:npoly-1].value = 1.D-6
      if (nback GT 0) then parinfo[npoly:npoly+nback-1].value = 1.D-6
   endelse
   ftol = 1d-20
   gtol = 1d-20
   xtol = 1d-20

;   ; Add constraints that the multiplicative term must be in the
;   ; bounds set by LIMITS, at least at the end points.
;   limits = [0.1, 10]
;   xmin = min(xvector)
;   xmax = max(xvector)
;   if (npoly EQ 1) then begin
;      parinfo[0].limits = limits
;   endif else if (npoly EQ 2) then begin
;      tiepar = replicate(parinfo1, 2)
;      tiepar[0].tied = 'P(0) + ' + string(xmin) + ' * P(1)'
;      tiepar[1].tied = 'P(0) + ' + string(xmax) + ' * P(1)'
;      tiepar[0].limits = limits
;      tiepar[1].limits = limits
;      parinfo = [parinfo, tiepar]
;   endif else if (npoly EQ 3) then begin
;      tiepar[0].tied = 'P(0) + ' + string(xmin) + ' * P(1) + ' $
;       + string(xmin^2) + ' * P(2) * P(2)'
;      tiepar[1].tied = 'P(0) + ' + string(xmax) + ' * P(1) + ' $
;       + string(xmax^2) + ' * P(2) * P(2)'
;      tiepar[0].limits = limits
;      tiepar[1].limits = limits
;      parinfo = [parinfo, tiepar]
;   endif

   acoeff = mpfit('solve_poly_chi_fn', parinfo=parinfo, perror=perror, $
    maxiter=maxiter, ftol=ftol, gtol=gtol, xtol=xtol, $
    niter=niter, status=status, /quiet)

   if (arg_present(totchi2)) then $
    totchi2 = total( (solve_poly_chi_fn(acoeff))^2 )

   yfit = solve_poly_fn(acoeff, ymult=ymult, yadd=yadd)

   return
end

;------------------------------------------------------------------------------
; Fit for BFLUX = polynomial * AFLUX + background.
; The errors are the simple quadrature sum from AIVAR and BIVAR.
pro solve_poly_ratio1, xvector, aflux, bflux, aivar, $
 npoly=npoly1, nback=nback1, yfit=yfit, ymult=ymult, yadd=yadd, acoeff=acoeff, $
 totchi2=totchi2

   if (keyword_set(npoly1)) then npoly = npoly1[0] $
    else npoly = 1
   if (n_elements(nback1) GT 0) then nback = nback1[0] $
    else nback = 0
   wconstrain = 10 ; The weight of the constraints

   ; Construct the polynomial and additive vectors
   aarr = solve_poly_vectors(xvector, npoly)
   barr = solve_poly_vectors(xvector, nback)

   npix = n_elements(xvector)
   nterm = npoly + nback

   ; Compute the errors to use
   ; If the problem is ill-constrained, then try to make it constrained
   ; by adding 1 to all the inverse sigmas.  Note that this will still
   ; be underconstrained if there are too few data points, but the SVD
   ; algorithm should still do something sensible.
   sqivar = sqrt(aivar)
   if (total(sqivar GT 0) LT nterm) then sqivar = sqivar + 1.

   ; Construct the empty matrices
   mmatrix = dblarr(npix+nterm, nterm)
   bvec = dblarr(npix+nterm)

   ; Put the measured values into the matrix
   for i=0, npoly-1 do $
    mmatrix[0:npix-1,i] = aflux[*] * sqivar * aarr[*,i]
   for j=0, nback-1 do $
    mmatrix[0:npix-1,npoly+j] = barr[*,j]
   bvec[0:npix-1] = bflux[*] * sqivar

   ; Put the constraints into the matrix: a[0] = 1
   mmatrix[npix,0] = 1 * wconstrain
   bvec[npix] = 1 * wconstrain

   ; Put the constraints into the matrix: a[1...] = 0
   if (npoly GT 1) then begin
      for i=1, npoly-1 do mmatrix[npix+i,i] = 1 * wconstrain
      bvec[npix+1:npix+npoly-1] = 0 * wconstrain
   endif

   ; Put the constraints into the matrix: b[0...] = 0
   if (nback GT 0) then begin
      for j=0, nback-1 do mmatrix[npix+npoly+j,npoly+j] = 1 * wconstrain
   endif

   ; Now invert the matrix
   mmatrixt = transpose(mmatrix)
   mm = mmatrixt # mmatrix
   if (nterm EQ 1) then begin
      mmi = 1.0 / (mm + (mm EQ 0))
   endif else begin
      svdc, mm, ww, uu, vv, /double
      mmi = 0 * vv
      ; The last term below is to protect against divide-by-zero
      ; in the degenerate case.
      for i=0L, nterm-1 do mmi[i,*] = vv[i,*] / (ww[i] + (ww[i] EQ 0))
      mmi = mmi ## transpose(uu)
   endelse

   acoeff = mmi # (mmatrixt # bvec)
   totchi2 = total( (mmatrix # acoeff - bvec)^2, /double )

   ymult = acoeff[0] * aarr[*,0]
   for i=1, npoly-1 do ymult = ymult + acoeff[i] * aarr[*,i]

   yadd = 0 * aflux[*]
   if (nback GT 0) then yadd = yadd + acoeff[npoly:npoly+nback-1] ## barr

   yfit = ymult * aflux[*] + yadd

   return
end
;------------------------------------------------------------------------------
pro solve_poly_ratio, xvector, aflux, bflux, aivar, bivar, $
 npoly=npoly, nback=nback, yfit=yfit, ymult=ymult, yadd=yadd, $
 acoeff=acoeff, totchi2=totchi2, $
 inparams=inparams1, status=status, perror=perror

   if (n_params() EQ 4) then begin
      solve_poly_ratio1, xvector, aflux, bflux, aivar, $
       npoly=npoly, nback=nback, yfit=yfit, ymult=ymult, yadd=yadd, $
       acoeff=acoeff, totchi2=totchi2
   endif else if (n_params() EQ 5) then begin
      if (keyword_set(inparams1)) then begin
         inparams = inparams1
      endif else begin
         ; Do a 1st pass re-fitting without scaling the errors
         ; in order to get the initial guess for the nonlinear fitting
         thisivar = 0 * aivar
         indx = where(aivar GT 0 AND bivar GT 0, ct)
         if (ct GT 0) then $
          thisivar[indx] = 1. / (1./aivar[indx] + 1./bivar[indx])
         solve_poly_ratio, xvector, aflux, bflux, thisivar, $
          npoly=npoly, nback=nback, acoeff=inparams
      endelse
      solve_poly_ratio2, xvector, aflux, bflux, aivar, bivar, $
       npoly=npoly, nback=nback, yfit=yfit, ymult=ymult, yadd=yadd, $
       acoeff=acoeff, totchi2=totchi2, $
       inparams=inparams, status=status, perror=perror
   endif else begin
      message, 'Invalid inputs!'
   endelse

   return
end
;------------------------------------------------------------------------------
