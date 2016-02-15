;------------------------------------------------------------------------------
;+
; NAME:
;   ext_ccm
;
; PURPOSE:
;   Return extinction curve from CCM (1989), defined in the wavelength
;   range [1250,33333] Angstroms.
;
; CALLING SEQUENCE:
;   Alam = ext_ccm( lambda, [ Rv ] )
;
; INPUTS:
;   lambda:      Wavelength(s) in Angstroms
;
; OPTIONAL INPUTS:
;   Rv:          Value of R_V; default to 3.1
;
; OUTPUTS:
;   Alam:        Return value A(lambda)/A(V)
;
; COMMENTS:
;   Reference is Cardelli, J.A., Clayton, G.C., & Mathis, J.S. 1989,
;   ApJ, 345, 245-256.
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Written D. Schlegel, 8 September 1997, Durham
;-
;------------------------------------------------------------------------------
function ext_ccm, lambda, Rv

   if (n_params() LT 1) then begin
      message, 'Syntax - Alam = ext_ccm( lambda, [ Rv ] )'
   endif

   if (NOT keyword_set(Rv)) then Rv = 3.1

   xx = 10000.d0 / lambda
   Alam = fltarr(n_elements(lambda))

   ; Limits for CCM fitting function
   qLO = where(xx GT 8.0)                ; No data, lambda < 1250 Ang
   qUV = where(xx GT 3.3 AND xx LE 8.0)  ; UV + FUV
   qOPT = where(xx GE 1.1 AND xx LE 3.3) ; Optical/NIR
   qIR = where(xx GT 0.3 AND xx LT 1.1)  ; IR
   qHI = where(xx LT 0.3)                ; No data, lambda > 33,333 Ang

   ; For lambda < 1250 Ang, arbitrarily return Alam=5
   if (qLO[0] NE -1) then begin
      Alam[qLO] = 5.0
   endif

   if (qUV[0] NE -1) then begin
      xt = xx[qUV]
      afac = 1.752 - 0.316*xt - 0.104 / ( (xt-4.67)^2 + 0.341 )
      bfac = -3.090 + 1.825*xt + 1.206 / ( (xt-4.62)^2 + 0.263 )
      Fa = -0.04473*(xt-5.9)^2 - 0.009779*(xt-5.9)^3
      Fb = 0.2130*(xt-5.9)^2 + 0.1207*(xt-5.9)^3
      qq = where(xt GE 5.9 AND xt LE 8.0)
      if (qq[0] NE -1) then begin
         afac[qq] = afac[qq] + Fa
         bfac[qq] = bfac[qq] + Fb
      endif
      Alam[qUV] = afac + bfac / Rv
   endif

   if (qOPT[0] NE -1) then begin
      yy = xx[qOPT] - 1.82
      afac = 1.0 + 0.17699*yy - 0.50447*yy^2 - 0.02427*yy^3 + 0.72085*yy^4 $
       + 0.01979*yy^5 - 0.77530*yy^6 + 0.32999*yy^7
      bfac = 1.41338*yy + 2.28305*yy^2 + 1.07233*yy^3 - 5.38434*yy^4 $
       - 0.62251*yy^5 + 5.30260*yy^6 - 2.09002*yy^7
      Alam[qOPT] = afac + bfac / Rv
   endif

   if (qIR[0] NE -1) then begin
      yy = xx[qIR]^1.61
      afac = 0.574*yy
      bfac = -0.527*yy
      Alam[qIR] = afac + bfac / Rv
   endif

   ; For lambda > 33,333 Ang, arbitrarily extrapolate the IR curve
   if (qHI[0] NE -1) then begin
      yy = xx[qHI]^1.61
      afac = 0.574*yy
      bfac = -0.527*yy
      Alam[qHI] = afac + bfac / Rv
   endif

   return, Alam
end
;------------------------------------------------------------------------------
