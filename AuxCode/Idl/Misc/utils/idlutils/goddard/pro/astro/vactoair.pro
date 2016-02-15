pro vactoair,wave
;+
; NAME:
;	VACTOAIR
; PURPOSE:
;	Convert vacuum wavelengths to air wavelengths
; EXPLANATION:
;	Corrects for the index of refraction of air under standard conditions.  
;	Wavelength values below 2000 A will not be altered.  Accurate to 
;	about 0.005 A 
;
; CALLING SEQUENCE:
;	VACTOAIR, WAVE
;
; INPUT/OUTPUT:
;	WAVE - Wavelength in Angstroms, scalar or vector
;		WAVE should be input as vacuum wavelength(s), it will be
;		returned as air wavelength(s).  WAVE is always converted to
;		double precision
;
; EXAMPLE:
;	If the vacuum wavelength is  W = 2000, then 
;
;	IDL> VACTOAIR, W 
;
;	yields an air wavelength of W = 1999.353 Angstroms
;
; METHOD:
;	An approximation to the 4th power of inverse wavenumber is used
;	See IUE Image Processing Manual   Page 6-15.
;
; REVISION HISTORY
;	Written, D. Lindler 1982 
;	Documentation W. Landsman  Feb. 1989
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
  On_error,2
  if N_params() EQ 0 then begin
     print,'Syntax - VACTOAIR, Wave'
     return
  endif

  wave2 = double(wave)*wave 
  fact = 1.0 + 2.735182e-4 + 131.4182/wave2 + 2.76249e8/(wave2*wave2)
  fact = fact * ( wave GE 2000. ) + 1.0*( wave LT 2000.0 )

; Convert wavelengths

  wave = wave/fact

  return
  end                        
