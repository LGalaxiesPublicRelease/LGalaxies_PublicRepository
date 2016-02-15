;------------------------------------------------------------------------------
;+
; NAME:
;   ext_odonnell
;
; PURPOSE:
;   Return extinction curve from Odonnell (1994), defined in the wavelength
;   range [3030,9091] Angstroms.  Outside this range, use CCM (1989).
;
; CALLING SEQUENCE:
;   Alam = ext_odonnell( lambda, [ Rv ] )
;
; INPUTS:
;   lambda:      Wavelength(s) in Angstroms
;
; OPTIONAL INPUTS:
;   Rv:          Value of R_V (scalar); default to 3.1
;
; OUTPUTS:
;   Alam:        Return value A(lambda)/A(V)
;
; COMMENTS:
;   Reference is O'Donnell, J. E. 1994, ApJ, 422, 158-163.
;
; PROCEDURES CALLED:
;   ext_ccm()
;
; REVISION HISTORY:
;   Written D. Schlegel, 8 September 1997, Durham
;-
;------------------------------------------------------------------------------
function ext_odonnell, lambda, Rv

   if (n_params() LT 1) then begin
      message, 'Syntax - Alam = ext_odonnell( lambda, [ Rv ] )'
   endif

   if (NOT keyword_set(Rv)) then Rv = 3.1

   xx = 10000.d0 / lambda
   Alam = fltarr(n_elements(lambda))
 
   ; Limits for CCM fitting function
   qOPT = where(xx GE 1.1 AND xx LE 3.3) ; Optical/NIR fitting function
   qBAD = where(xx LT 1.1 OR xx GT 3.3)  ; Outside fitting range

   if (qOPT[0] NE -1) then begin
      yy = xx[qOPT] - 1.82
      afac = 1.0 + 0.104*yy - 0.609*yy^2 + 0.701*yy^3 + 1.137*yy^4 $
       - 1.718*yy^5 - 0.827*yy^6 + 1.647*yy^7 - 0.505*yy^8
      bfac = 1.952*yy + 2.908*yy^2 - 3.989*yy^3 - 7.985*yy^4 $
       + 11.102*yy^5 + 5.491*yy^6 - 10.805*yy^7 + 3.347*yy^8
      Alam[qOPT] = afac + bfac / Rv
   endif
 
   ; Outside the O'Donnell fitting range of [3030,9091] Angstroms,
   ; use the CCM fitting function.
   if (qBAD[0] NE -1) then begin
      Alam[qBAD] = ext_ccm(lambda[qBAD], Rv)
   endif

   return, Alam
end
;------------------------------------------------------------------------------
