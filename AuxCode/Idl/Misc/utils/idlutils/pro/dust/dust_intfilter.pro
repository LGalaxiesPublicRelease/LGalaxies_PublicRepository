;------------------------------------------------------------------------------
;+
; NAME:
;   dust_intfilter
;
; PURPOSE:
;   Integrate the extinction curve over a given source function + filter curve.
;
; CALLING SEQUENCE:
;   avec = dust_intfilter( ixval, ffilename=, sfilename=, [ atmfilename=, $
;    rv=, anorm=, dlam=, zsource=, /debug  ] )
;
; INPUTS:
;   ixval      - Temperature-corrected 100-micron flux (IX) from SFD maps,
;                in MJy/sr (scalar or vector).
;
; REQUIRED KEYWORDS:
;   ffilename  - ASCII file with wavelengths in Ang (1st column) and throughput
;                (2nd column) for the filter response curve.
;   sfilename  - ASCII file with wavelengths in Ang (1st column) and flux in
;                f_lambda (2nd column) for the source function.  Note that
;                since we assume f_lambda (flux/Ang), we multiply this
;                by one power of the wavelength to convert to a flux in
;                photon number per Angstrom.
;
; OPTIONAL KEYWORDS:
;   atmfilename- ASCII file with wavelengths in Ang (1st column) and magnitudes
;                of atmospheric extinction (2nd column).  If specified,
;                then the filter response is multiplied by 10^(-ATM/2.5).
;   rv         - Extinction curve parameter R_V; default to 3.1.
;   anorm      - Normalization factor for multiplying IX to obtain
;                extinction at 1 micron; default to SFD normalization
;                of (1.319)*(0.0184) mag/(MJy/sr).
;   dlam       - Spacing of numeric integration in Angstroms; default
;                to 1 Ang.
;   zsource    - Redshift of source; default to 0.
;   debug      - If set, then make debugging plots with SPLOT.
;
; OUTPUTS:
;   avec       - Extinction in magnitudes (same dimensions as IXVAL).
;
; COMMENTS:
;   The Galactic extinction curve is that from O'Donnell (1994)
;   and Cardelli, Clayton & Mathis (1989).
;
;   The integrations (and optional debugging plots) are limited to
;   the wavelength range within which the filter curve is positive-valued.
;
;   The filter curve, source function, and atmospheric extinction curve
;   are cached between calls.  If the same file is specified on a
;   subsequent call, then those cached values are used.
;
; PROCEDURES CALLED:
;   ext_odonnell()
;   numlines()
;   readcol
;
; DATA FILES:
;
; REVISION HISTORY:
;   01-Dec-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
function dust_intfilter, ixval, ffilename=ffilename1, sfilename=sfilename1, $
 atmfilename=atmfilename1, rv=rv, anorm=anorm, dlam=dlam1, $
 zsource=zsource1, debug=debug

   common com_dust_intfilter, ffilename, fwave, fthru, $
    sfilename, swave, sflux, zsource, atmfilename, awave, aext, $
    dlam, lambda, f_tabulate, s_tabulate, atm_tabulate

   if (n_params() LT 1) then $
    message, 'Must specify IXVAL'
   if (NOT keyword_set(ffilename1)) then $
    message, 'FFILENAME must be specified'
   if (NOT keyword_set(sfilename1)) then $
    message, 'SFILENAME must be specified'
   if (NOT keyword_set(rv)) then rv = 3.1
   if (NOT keyword_set(anorm)) then anorm = 1.319 * 0.0184

   if (NOT keyword_set(dlam)) then dlam = 0
   if (NOT keyword_set(dlam1)) then dlam1 = 1. ; default value
   if (dlam1 NE dlam) then begin
      dlam = dlam1
      qnewlam = 1
   endif

   if (NOT keyword_set(zsource)) then zsource = 0
   if (keyword_set(zsource1)) then begin
      if (zsource1 NE zsource) then qnewz = 1
      zsource = zsource1
   endif else begin
      if (zsource NE 0) then qnewz = 1
      zsource = 0
   endelse

   ;----------
   ; Read in the filter curve

   if (NOT keyword_set(ffilename)) then ffilename = ''
   if (ffilename1 NE ffilename) then begin
      ffilename = ffilename1
      if (numlines(ffilename) LE 0) then begin
         print, 'Filter curve file is empty ' + ffilename
         ffilename = ''
         return, 0
      endif
      readcol, ffilename, fwave, fthru, /silent, format='D,D'
      qnewffile = 1
      qnewlam = 1
   endif
   if (NOT keyword_set(fwave) OR NOT keyword_set(fthru)) then $
    message, 'Filter response function not set'

   ;----------
   ; Read in the source function

   if (NOT keyword_set(sfilename)) then sfilename = ''
   if (sfilename1 NE sfilename) then begin
      sfilename = sfilename1
      if (numlines(sfilename) LE 0) then begin
         print, 'Source function file is empty ' + sfilename
         sfilename = ''
         return, 0
      endif
      readcol, sfilename, swave, sflux, /silent, format='D,D'
      qnewsfile = 1
   endif
   if (NOT keyword_set(swave) OR NOT keyword_set(sflux)) then $
    message, 'Source function not set'

   ;----------
   ; Read in the optional atmospheric extinction file

   if (NOT keyword_set(atmfilename)) then atmfilename = ''
   if (keyword_set(atmfilename1)) then begin
      if (atmfilename1 NE atmfilename) then begin
         atmfilename = atmfilename1
         if (numlines(atmfilename) LE 0) then begin
            print, 'Atmospheric extinction file is empty ' + atmfilename
            atmfilename = ''
            awave = 0
            aext = 0
            return, 0
         endif
         readcol, atmfilename, awave, aext, /silent, format='D,D'
         qnewafile = 1
      endif
   endif else begin
      atmfilename = ''
      atm_tabulate = 1 ; Default to no atmospheric extinction
   endelse

   ;----------
   ; Select a set of wavelengths for the integration.

   if (keyword_set(qnewlam) OR keyword_set(qnewffile)) then begin
      indx = where(fthru GT 0, ct)
      if (ct EQ 0) then $
       message, 'No positive values in filter response curve'
      wavemin = floor( min(fwave[indx]) )
      wavemax = max(fwave[indx])
      nwave = ceil((wavemax - wavemin) / dlam ) + 1
      lambda = wavemin + dindgen(nwave) * dlam
   endif

   ;----------
   ; Evaluate the filter curve at the integration points

   if (keyword_set(qnewlam) OR keyword_set(qnewffile) $
    OR NOT keyword_set(f_tabulate)) then begin
      linterp, fwave, fthru, lambda, f_tabulate, missing=0.0
      f_tabulate = f_tabulate / max(f_tabulate) ; Arbitrary normalization
   endif

   ;----------
   ; Evaluate the source function at the integration points.
   ; Multiply by one power of wavelength to convert from f_lambda
   ; (i.e., erg/s/cm^2/Ang) to photon number flux (i.e., photons/Ang).
   ; In the photon-number units, the R-J tail of a blackbody falls off as
   ; lambda^(-1).

   if (keyword_set(qnewlam) OR keyword_set(qnewsfile) OR keyword_set(qnewz) $
    OR NOT keyword_set(s_tabulate)) then begin
      linterp, swave*(1.+zsource), sflux, lambda, s_tabulate
      s_tabulate = s_tabulate * (lambda / lambda[0])
      s_tabulate = s_tabulate / max(s_tabulate) ; Arbitrary normalization
   endif

   ;----------
   ; Evaluate the atmospheric extinction at the integration points

   if (keyword_set(qnewafile) $
    OR (keyword_set(atmfilename) AND NOT keyword_set(atm_tabulate)) $
    OR (keyword_set(qnewlam) AND keyword_set(atmfilename)) ) then begin $
      linterp, awave, aext, lambda, atm_tabulate
      atm_tabulate = 10.d0^(-0.4 * atm_tabulate)
   endif

   ;----------
   ; Perform the integral for each value of IXVAL.

   amag = 0. * ixval ; Allocate memory for return value
   fsvec = f_tabulate * s_tabulate * atm_tabulate
   denom = total(fsvec)
   aratio = ext_odonnell(lambda, rv) / (ext_odonnell(10000., rv))[0]
   for i=0L, n_elements(ixval)-1 do begin
      avec = anorm * ixval[i] * aratio
      numer = total(fsvec * 10^(-avec/2.5))
      amag[i] = -2.5 * alog10(numer / denom)

      if (keyword_set(debug)) then begin
         csize = 1.5
         splot, lambda, f_tabulate, color='blue', $
          xtitle='Wavelength [Ang]', charsize=csize
         soplot, lambda, s_tabulate, color='green'
         soplot, lambda, 10^(-avec/2.5), color='red'
         if (n_elements(atm_tabulate) GT 1) then $
          soplot, lambda, atm_tabulate, color='magenta'
         soplot, lambda, fsvec * 10^(-avec/2.5)
         sxyouts, !x.crange[0], 0.90, charsize=csize, $
          '  Filter response', color='blue'
         sxyouts, !x.crange[0], 0.85, charsize=csize, $
          '  Source fn (photons/Ang)', color='green'
         sxyouts, !x.crange[0], 0.80, charsize=csize, $
          '  Extinction', color='red'
         if (n_elements(atm_tabulate) GT 1) then $
          sxyouts, !x.crange[0], 0.75, charsize=csize, $
           '  Atmosphere', color='magenta'
         sxyouts, !x.crange[0], 0.70, charsize=csize, $
          '  Filter * Source * Extinction' $
          + (n_elements(atm_tabulate) GT 1 ? ' * Atmosphere' : '')
         print, 'Press any key...'
         cc = get_kbrd(1)
      endif
   endfor

   return, amag
end
;------------------------------------------------------------------------------
