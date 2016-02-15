pro airtovac,wave                   
;+
; NAME:
;       AIRTOVAC
; PURPOSE:
;       Convert air wavelengths to vacuum wavelengths 
; EXPLANATION:
;       Wavelengths are corrected for the index of refraction of air under 
;       standard conditions.  Wavelength values below 2000 A will not be 
;       altered.  Uses the IAU standard for conversion given in Morton 
;       (1991 Ap.J. Suppl. 77, 119)
;
; CALLING SEQUENCE:
;       AIRTOVAC, WAVE
;
; INPUT/OUTPUT:
;       WAVE - Wavelength in Angstroms, scalar or vector
;               WAVE should be input as air wavelength(s), it will be
;               returned as vacuum wavelength(s).  WAVE is always converted to
;               double precision upon return.
;
; EXAMPLE:
;       If the air wavelength is  W = 6056.125 (a Krypton line), then 
;       AIRTOVAC, W yields an vacuum wavelength of W = 6057.8019
;
; METHOD:
;       See Morton (Ap. J. Suppl. 77, 119) for the formula used
;
; REVISION HISTORY
;       Written W. Landsman                November 1991
;       Converted to IDL V5.0   W. Landsman   September 1997
;-
   On_error,2

  if N_params() EQ 0 then begin
      print,'Syntax - AIRTOVAC, WAVE'
      print,'WAVE (Input) is the air wavelength in Angstroms'
      print,'On output WAVE contains the vacuum wavelength in Angstroms'
      return
  endif

  sigma2 = (1d4/double(wave) )^2.              ;Convert to wavenumber squared

; Compute conversion factor

  fact = 1.D + 6.4328D-5 + 2.94981D-2/(146.D0 - sigma2) + $
                            2.5540D-4/( 41.D0 - sigma2)
    
  fact = fact*(wave GE 2000.) + 1.0*(wave LT 2000.0)

  wave = wave*fact              ;Convert Wavelength

  return            
  end
