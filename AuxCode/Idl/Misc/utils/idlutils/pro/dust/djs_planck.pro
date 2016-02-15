function djs_planck, Ttemp, nu_or_lambda, dBdT, units=units, mjy=mjy, wcm2=wcm2
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     DJS_PLANCK returns the spectral radiance of a blackbody.
;
;DESCRIPTION:  
;    IDL function to return the spectral radiance of a blackbody,
;    i.e. the Planck curve, in units of either MJy/steradian (I_nu)
;    or watts/cm^2/steradian (nu_I_nu).
;    The blackbody temperature and either frequency (in icm or GHz)
;    or wavelength (in microns) are inputs to the function.  The
;    routine also optionally returns the derivative with respect to 
;    temperature, in units of MJy/sr/K or W/cm^2/sr/K.
;
;CALLING SEQUENCE:  
;     RESULT = DJS_PLANCK (temperature, nu_or_lambda [,dBdT] $
;              [,UNITS=units], [/MJY], [/WCM2])
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     RESULT        O   flt [arr]  Spectral radiance at each wavelength. 
;                                  Units: W/cm^2/sr/K if /WCM2 specified
;                                         MJy/sr      if /MJY specfied
;     TEMPERATURE   I   flt        Temperature of blackbody, in K.
;     NU_OR_LAMBDA  I   flt        Frequency or wavelength at which to 
;                                  calculate spectrum. Units are as 
;                                  specified with UNITS keyword.
;     dBdT         [O]  flt [arr]  Derivative of Planck with respect to 
;                                  temperature. 
;     UNITS        [I]  str        'Microns', 'icm', or 'GHz' to 
;                                  identify units of NU_OR_LAMBDA. Only 
;                                  first character is required.  If 
;                                  left out, default is 'microns'.
;     /MJY          I   key        Sets output units to MJy/sr
;     /WCM2         I   key        Sets output units to W/cm^2/sr
;
;WARNINGS:
;     1.  One of /MJY or /WCM2 MUST be specified.  
;     2.  Routine gives incorrect results for T < 1 microKelvin and
;            wavelengths shortward of 1.e-10 microns.  (So sue me).
;
;EXAMPLE:
;     To produce a 35 K spectrum in MJy/sr at 2, 4, 6, 8, 10 microns:
;
;       wavelength = 2. + 2.*findgen(5)
;       temp = 35.
;       blackbody = djs_planck(temp, wavelength, units='micron', /mjy)
;
;     One could also get back the derivative by including it in the
;     call:
;       blackbody = djs_planck(temp, wavelength, deriv, units='m', /mjy)
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;     Identifies units using the UNITS keyword, then converts the 
;     supplied independent variable into microns to evaluate the 
;     Planck function.  Uses Rayleigh-Jeans and Wien approximations 
;     for the low- and high-frequency end, respectively.  Reference: 
;     Allen, Astrophysical Quantities, for the Planck formula.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     None
;  
;MODIFICATION HISTORY:
;    Written by Rich Isaacman, General Sciences Corp.  17 April 1991
;    Revised by RBI 27 Jan 1992 to use updated fundamental constants 
;         (SPR 9449)
;    Revised by RBI 29 Jan 1992 to calculate derivatives only when 
;         necessary
;    Revised by RBI 3 Feb 1992 to redimension output to a scalar if only 
;       a single wavelength is supplied  (SPR 9459)
;    Revised by RBI 6 Mar 92 to return either MJy/sr or (!) W/cm^2/sr
;    Revised by RBI 1 Jun 92 to fix single-wavelength failure when no
;       derivative is requested (SPR 9738), and to use MESSAGE.
;    RBI corrected error in derivative calculation SPR 9817 (17 Jul 92)
;    RBI corrected error in Wien and RJ tails SPR 10392 (24 Dec 92)
;	 but didn't get it quite right (Piper/Kryszak, 28-Dec-92)
;    Revised by David Schlegel 10-Mar-1999 to allow calling with temperature
;        and/or wavelength as a vector; converted to IDL-5; renamed from
;        PLANCK() to DJS_PLANCK().  Note that this code was copied from
;        the COBE analysis software.
;
; SPR 9616
;.TITLE
; Routine DJS_PLANCK
;-
;
; Check on input parameters
;
on_error, 2
if n_elements(nu_or_lambda) lt 1 or n_elements(Ttemp) lt 1 then $
     message,'CALLING SEQUENCE: bbflux = djs_planck (temp,wavelength,units=<units>)'
if not keyword_set(mjy) and not keyword_set(wcm2) then $
     message, 'Either /MJy or /Wcm2 must be specified!'
if keyword_set(mjy) and keyword_set(wcm2) then $
     message, 'Only one of /MJy or /Wcm2 may be specified!'
;
makederiv = n_params() eq 3           ; see whether dBdT is requested
if n_elements(units) lt 1 then begin
   units = 'm'
;   message, /continue, "Wavelength units are assumed to be microns."
endif
T = float(ttemp) > 1.e-06                ;force temperature to be real*4
;
; Define some necessary constants
;
c = 299792.458d0                         ;speed of light, Physics Today Aug 90
hck = 14387.69d0                         ;h*c/k             "      "
thcc = 1.1910439d0                       ;2*h/c^2           "     "  
coeff = thcc/hck^4 * 1.e04               ;Stephan-Boltzmann * 15/pi^5
;
; Convert nu_or_lambda into lambda (in microns) depending on specified units
;
units0 = strupcase (strmid (units,0,1))
case units0 of
   'M':   lambda = nu_or_lambda > 1.e-10                ; microns
   'I':   lambda = 1.e04 / (nu_or_lambda > 1.e-10)      ; icm, avoiding nu=0.
   'G':   lambda = c / (nu_or_lambda > 1.e-10)          ; GHz, avoiding nu=0.
   else:  message, "Units must be specified as 'microns', 'icm', or 'GHz'"
endcase
; 
;  Variable fhz is a scale factor used to go from units of nu_I_nu units to 
;  MJy/sr if the keyword is set.  In that case, its value is
;  frequency in Hz * w/cm^2/Hz ==> MJy
;
if keyword_set(wcm2) then fhz = 1. + fltarr(n_elements(lambda))
if keyword_set(mjy) then fhz = c/lambda * 1.e-15         
;
;  Introduce dimensionless variable chi, used to check whether we are on 
;  Wien or Rayleigh Jeans tails
;
chi = hck / lambda / T
val = fltarr(n_elements(chi))
if makederiv then dBdT = fltarr(n_elements(chi))
;
;  Start on Rayleigh Jeans side
;
rj = where (chi lt 0.001)
if rj[0] ne -1 then begin
    val[rj] = coeff * T^4 * chi[rj]^3 / fhz[rj]
    if makederiv then dBdT[rj] = val[rj] / T
endif
;
;  Now do nonapproximate part
;
exact = where (chi ge 0.001 and chi le 50)
if exact[0] ne -1 then begin
    chi_ex = chi[exact]
    val[exact] = coeff * T^4 * chi_ex^4 / (exp(chi_ex) - 1.) / fhz[exact]
    if makederiv then dBdT[exact] = $
	val[exact] * chi_ex / T / (1. - exp(-chi_ex))
endif
;
;  ...and finally the Wien tail
;
wien = where (chi gt 50.)
if wien[0] ne -1 then begin
    chi_wn = chi[wien]
    val[wien] = coeff * T^4 * chi_wn^4 * exp(-chi_wn) / fhz[wien]
    if makederiv then dBdT[wien] = $
	val[wien] * chi_wn / T / (1. - exp(-chi_wn))
endif
;
;  Redimension to a scalar if only 1 temperature and 1 wavelength supplied
;  Modified by DJS 10-Mar-1999
;
if n_elements(val) eq 1 then begin
    if makederiv then dBdT = dBdT[0]
    val = val[0]
endif
return, val
;
end
;DISCLAIMER:
;
;This software was written at the Cosmology Data Analysis Center in
;support of the Cosmic Background Explorer (COBE) Project under NASA
;contract number NAS5-30750.
;
;This software may be used, copied, modified or redistributed so long
;as it is not sold and this disclaimer is distributed along with the
;software.  If you modify the software please indicate your
;modifications in a prominent place in the source code.  
;
;All routines are provided "as is" without any express or implied
;warranties whatsoever.  All routines are distributed without guarantee
;of support.  If errors are found in this code it is requested that you
;contact us by sending email to the address below to report the errors
;but we make no claims regarding timely fixes.  This software has been 
;used for analysis of COBE data but has not been validated and has not 
;been used to create validated data sets of any type.
;
;Please send bug reports to CGIS@COBECL.DNET.NASA.GOV.


