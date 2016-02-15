; INPUT
;  astr  : basic astrometry header (should be RA, DEC TAN)
; OPTIONAL INPUT
;  quad  : quadratic radial coefficient
; OUTPUT:
;  gsa   : GSSS astrometry structure. 


PRO astr2gsss, astr, gsa, cubecoeff=cubecoeff, quintcoeff=quintcoeff

  IF keyword_set(cubecoeff)  EQ 0 THEN cubecoeff  = 0.d
  IF keyword_set(quintcoeff) EQ 0 THEN quintcoeff = 0.d
  
  gsa = {gsss_astrometry, CTYPE: strarr(2), XLL:0, YLL:0, XSZ:0.D, YSZ:0.D, $
         PPO3:0.0D, PPO6:0.0D, CRVAL: dblarr(2), PLTSCL:0.0D, $
         AMDX:dblarr(13), AMDY:dblarr(13) }
  
  gsa.xll = 0
  gsa.yll = 0

  gsa.crval = astr.crval
  gsa.ctype = ['RA---GSS', 'DEC--GSS']

  gsa.xsz =  -1000.d
  gsa.ysz =  1000.d
  gsa.ppo3 = gsa.xsz*(astr.crpix[0]-0.5)  ; treat orientation coeffs as CRPIX
  gsa.ppo6 = gsa.ysz*(astr.crpix[1]-0.5)

  gsa.amdx = dblarr(13)
  gsa.amdy = dblarr(13)

; Map CD matrix onto coefficients, including cubic and quintic terms

  cd = astr.cd
  gsa.amdx[0]  = cd[0, 0] *3600. ; x  coeff
  gsa.amdx[1]  = cd[0, 1] *3600. ; y  
  gsa.amdx[8]  = cubecoeff * gsa.amdx[1] ; y x^2
  gsa.amdx[10] = cubecoeff * gsa.amdx[1] ; y y^2
  gsa.amdx[11] = cubecoeff * gsa.amdx[0] ; x r^2

  gsa.amdy[0]  = cd[1, 1] *3600. ; y  coeff (not a typo!)
  gsa.amdy[1]  = cd[1, 0] *3600. ; x
  gsa.amdy[8]  = cubecoeff * gsa.amdy[1] ; x y^2
  gsa.amdy[10] = cubecoeff * gsa.amdy[1] ; x x^2
  gsa.amdy[11] = cubecoeff * gsa.amdy[0] ; y r^2

  gsa.pltscl = 1.d

  return
END

