;------------------------------------------------------------------------------
;+
; NAME:
;   predict_thermal
;
; PURPOSE:
;   Read predictions of thermal dust emission from Finkbeiner et al. maps
;   and return sky intensity.
;
; CALLING SEQUENCE:
;   value = predict_thermal( [ gall, galb, nu=nu, infile=infile, $
;    skipline=skipline, outfile=outfile, resolution=resolution, model=model, $
;    interp=interp, noloop=noloop, verbose=verbose, ipath=ipath, units=units ] )
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   gall:       Galactic longitude(s) in degrees
;   galb:       Galactic latitude(s) in degrees
;   nu:         Frequency in GHz.  If this is a vector, it must be the same
;               dimension as GALL and GALB.  If this is a scalar, then it
;               applies to all positions GALL, GALB.
;   resolution: Set to one of the following (default is 'I4096'):
;               'I4096' : IRAS 4096^2 map (highest-resolution; default)
;               'I2048' : IRAS 2048^2 map
;               'I1024' : IRAS 1024^2 map
;               'D1024' : DIRBE 1024^2 map
;   model:      Model number (default to 8):
;               1: One-component, nu^1.5 emissivity
;               2: One-component, nu^1.7 emissivity
;               3: One-component, nu^2.0 emissivity
;               4: One-component, nu^2.2 emissivity
;               5: Two-component, alpha1=1.5, alpha2=2.6, Pollack et al. model
;               6: Two-component, both nu^2 emissivities, fit f+q
;               7: Two-component, alpha1=1.5, alpha2=2.6, fit f+q
;               8: Two-component, alpha1=1.67, alpha2=2.70, fit alphas+f+q
;   infile:     If set, then read GALL, GALB from this file.  If NU is not
;               set, then NU is read as the 3rd column of this same file.
;   skipline:   Number of lines to skip at the top of the input file
;   outfile:    If set, then write results to this file
;   interp:     Set this flag to return a linearly interpolated value
;               from the 4 nearest pixels.
;   noloop:     Set this flag to read all values at once without a FOR loop.
;               This is a faster option for reading a large number of values,
;               but requires reading an entire FITS image into memory.
;               (Actually, the smallest possible sub-image is read.)
;   verbose:    Set this flag for verbose output, printing pixel coordinates
;               and map values.  Setting NOLOOP disables this option.
;   ipath:      Path name for dust maps; default to path set by the
;               environment variable $DUST_DIR/map, or to the current
;               directory.
;   units:      Units for output values:
;               'MJy'    : MJy/sr (default)
;               'microK' : brightness (antenna) temperature [micro-Kelvin]
;               'thermo' : thermodynamic temperature [micro-Kelvin]
;                          assuming T(CMB) = 2.73 K
;
; OUTPUTS:
;   value:      Predicted emission value(s) from Lambert-projection maps
;
; COMMENTS:
;   These predictions are based upon the following paper:
;   "Extrapolation of Galactic Dust Emission at 100 Microns to CMBR
;   Frequencies using FIRAS", Finkbeiner, D. P., Davis, M., & Schlegel, D. J.,
;   ApJ, 1999, submitted (5 March 1999)
;
;   Either the coordinates GALL, GALB and NU must be set, or their values
;   must exist in the file INFILE.  Output is written to the variable VALUE
;   and/or the file OUTFILE.
;
; EXAMPLES:
;   Read the predicted thermal emission from dust at Galactic (l,b)=(12,+34.5) 
;   interpolating from the nearest 4 pixels, and return result in VALUE.
;   > value = predict_thermal(12, 34.5, nu=3000, /interp)
;
; PROCEDURES CALLED:
;   djs_planck()
;   fac_flux2temp()
;   planckcorr()
;   rdfloat
;   wcs_getval()
;
; DATA FILES:
;   FINK_Rmap_ngp.fits
;   FINK_Rmap_sgp.fits
;   SFD_d100_1024_ngp.fits
;   SFD_d100_1024_sgp.fits
;   SFD_i100_1024_ngp.fits
;   SFD_i100_1024_sgp.fits
;   SFD_i100_2048_ngp.fits
;   SFD_i100_2048_sgp.fits
;   SFD_i100_4096_ngp.fits
;   SFD_i100_4096_sgp.fits
;
; INTERNAL PROCEDURES:
;   kfactor()
;
; REVISION HISTORY:
;   10-Mar-1999  Written by David Schlegel, Princeton
;-
;------------------------------------------------------------------------------
; Return DIRBE color-correction (K-factor) for 100-micron map.
; "alpha" must be a scalar, but "temp" may be a vector.
function kfactor, alpha, temp

   ; We tabulate the K-factor for only the following emissivity profiles
   Kfitindx = [1.50, 1.67, 1.70, 2.00, 2.20, 2.60, 2.70]
   ia = (where(alpha EQ Kfitindx))[0]

   KfitAarr = $
    [ [  1.00000,  2.08243, -4.72422,  2.29118 ], $
      [  1.00000,  2.15146, -4.84539,  2.35210 ], $
      [  1.00000,  2.14106, -4.83639,  2.35919 ], $
      [  1.00000,  2.18053, -4.89849,  2.38060 ], $
      [  1.00000,  2.55941, -5.41290,  2.57867 ], $
      [  1.00000,  3.16383, -6.23131,  2.86900 ], $
      [  1.00000,  3.31600, -6.43306,  2.93939 ] ]

   KfitBarr = $
    [ [ -0.88339,  4.10104, -4.43324,  1.76240 ], $
      [ -0.87985,  4.10909, -4.43404,  1.76591 ], $
      [ -0.93625,  4.19278, -4.46069,  1.77103 ], $
      [ -0.80409,  3.95436, -4.27972,  1.70919 ], $
      [ -0.80318,  4.20361, -4.55598,  1.80207 ], $
      [ -0.50356,  4.07226, -4.70080,  1.87416 ], $
      [ -0.41568,  4.02002, -4.72432,  1.88865 ] ]
   nord = (size(KfitAarr))[1] - 1

; "temporary" function is used for better memory management
   log10T = alog10(temp)
   polya = KfitAarr[nord, ia]
   polyb = KfitBarr[nord, ia]
   for k=1, nord do begin
       polya = temporary(polya)*log10T + KfitAarr[nord-k, ia]
       polyb = temporary(polyb)*log10T + KfitBarr[nord-k, ia]
   endfor
log10T = 0 ; clear memory
   Kvalue = temporary(polya)/temporary(polyb)

   return, Kvalue
end
;------------------------------------------------------------------------------
function predict_thermal, gall, galb, nu=nu, infile=infile, $
 skipline=skipline, outfile=outfile, resolution=resolution, model=model, $
 interp=interp, noloop=noloop, verbose=verbose, ipath=ipath, units=units

   if (NOT keyword_set(infile) $
     AND (N_elements(gall) EQ 0 OR N_elements(galb) EQ 0 $
          OR N_elements(nu) EQ 0) $
    ) then begin
      print, 'Must set either coordinates and frequency or INFILE'
      return, -1
   endif
   if (NOT keyword_set(resolution)) then resolution = 'I4096'
   if (NOT keyword_set(model)) then model = 8
   if (NOT keyword_set(units)) then units = 'MJy'

   case resolution of
      'I4096': resname = 'i100_4096'
      'I2048': resname = 'i100_2048'
      'I1024': resname = 'i100_1024'
      'D1024': resname = 'd100_1024'
   endcase

   ; Set model parameters
   alpha1vec = [1.50, 1.70, 2.00, 2.20, 1.50, 2.00, 1.50, 1.67]
   alpha2vec = [0.00, 0.00, 0.00, 0.00, 2.60, 2.00, 2.60, 2.70]
   f1vec     = [1.00, 1.00, 1.00, 1.00, 0.25, 0.00261, 0.0309, 0.0363]
   q1q2vec   = [1.00, 1.00, 1.00, 1.00, 0.61, 2480.0, 11.2, 13.0]

   ; Rfita contains fit coefficients for T2_of_R.
   RfitAarr = $
    [[3.0998E+00, 2.7984E-01, 3.7562E-02, 7.5493E-03, 1.3661E-03, 1.1266E-04], $
     [3.0536E+00, 2.6701E-01, 3.3732E-02, 6.6142E-03, 1.2585E-03, 1.1106E-04], $
     [2.9897E+00, 2.4984E-01, 2.8390E-02, 5.3173E-03, 1.1161E-03, 1.1106E-04], $
     [2.9482E+00, 2.4061E-01, 2.5895E-02, 4.5086E-03, 9.9279E-04, 1.0790E-04], $
     [2.9206E+00, 2.3254E-01, 2.3506E-02, 4.0781E-03, 1.0048E-03, 1.2004E-04], $
     [2.9900E+00, 2.5041E-01, 2.9688E-02, 6.5641E-03, 1.5688E-03, 1.6542E-04], $
     [2.8874E+00, 2.4172E-01, 2.9369E-02, 4.7867E-03, 9.7237E-04, 1.1410E-04], $
     [2.8723E+00, 2.4071E-01, 2.9625E-02, 4.7196E-03, 9.3207E-04, 1.1099E-04] ]

; --- old incorrect array
;    [[2.9268E+00, 3.8419E-01, 5.0233E-02, 1.0852E-02, 3.0738E-03, 5.0595E-04], $
;     [2.8483E+00, 3.8044E-01, 4.6584E-02, 9.0938E-03, 2.7038E-03, 5.4664E-04], $
;     [2.7334E+00, 3.7537E-01, 4.1712E-02, 6.8839E-03, 2.0316E-03, 6.0311E-04], $
;     [2.6556E+00, 3.7377E-01, 3.9898E-02, 5.7662E-03, 1.4638E-03, 6.3723E-04], $
;     [2.9206E+00, 2.3254E-01, 2.3506E-02, 4.0781E-03, 1.0048E-03, 1.2004E-04], $
;     [2.9900E+00, 2.5041E-01, 2.9688E-02, 6.5641E-03, 1.5688E-03, 1.6542E-04], $
;     [2.8874E+00, 2.4172E-01, 2.9369E-02, 4.7867E-03, 9.7237E-04, 1.1410E-04], $
;     [2.8723E+00, 2.4071E-01, 2.9625E-02, 4.7196E-03, 9.3207E-04, 1.1099E-04] ]

   ; Zeta integrals for alpha=[1.50, 1.67, 1.70, 2.00, 2.20, 2.60, 2.70]
   ; from equn (15) of Finkbeiner et al.
   Zindx = [1.50, 1.67, 1.70, 2.00, 2.20, 2.60, 2.70]
   Zintegral = [5.3662E+01, 7.0562E+01, 7.4100E+01, 1.2208E+02, $
    1.7194E+02, 3.4855E+02, 4.1770E+02]

   ; Select parameters for this model
   alpha1 = alpha1vec[model-1]
   alpha2 = alpha2vec[model-1]
   f1 = f1vec[model-1]
   q1q2 = q1q2vec[model-1]
   nu100 = 2997.92458 ; Frequency in GHz for 100-microns
   RfitA = RfitAarr[*,model-1]

   ; If INFILE is defined, then read galactic coordinates from that file
   if (keyword_set(infile)) then begin
      if (keyword_set(nu)) then begin
         rdfloat, infile, gall, galb, skipline=skipline
      endif else begin
         rdfloat, infile, gall, galb, nu, skipline=skipline
         ; RDFLOAT leaves NU undefined if third column is missing. 
         if (N_elements(nu) EQ 0) then begin
            message, 'Must include NU in 3rd column or set as keyword'
            return, -1
         endif 
      endelse
   endif

   ; Convert nu to a float
   nu = float(nu)

   if keyword_set(ipath) then begin
       dum = findfile(ipath+'SFD*.fits', count=ct)
       if (ct EQ 0) then begin
           message, 'Bad file path!'
           return, -1
       endif
   endif

   if (NOT keyword_set(ipath)) then begin
      dust_dir = getenv('DUST_DIR')
      if (dust_dir NE '') then ipath = dust_dir+'/maps/' $
       else ipath = './'
   endif

   dum = findfile(ipath+'SFD*.fits', count=ct)
   if (ct EQ 0) then begin
       message, 'No data files found in path'
       return, -1
   endif

   ; Read the 100-micron map
   I100 = wcs_getval( $
    ['SFD_'+resname+'_ngp.fits', 'SFD_'+resname+'_sgp.fits'], $
    gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)
   lnR = alog( wcs_getval( $
    ['FINK_Rmap_ngp.fits', 'FINK_Rmap_sgp.fits'], $
    gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose) )

   if (model LE 4) then begin
      ; SINGLE-COMPONENT MODEL: Evaluate equn (1) from Finkbeiner et al

      ; Compute ln(T1) from ln(Rmap).  Loop through the coeffs to save memory.
      lnT1 = fltarr(N_elements(lnR))
      for ii=0, N_elements(RfitA)-1 do lnT1 = lnT1 + RfitA[ii] * lnR^ii
lnR = 0 ; clear memory
      T1 = exp(lnT1)
lnT1 = 0 ; clear memory

      Inu = I100 * (nu/nu100)^alpha1 * djs_planck(T1,nu,units='GHz',/MJy) / $
       (djs_planck(T1,nu100,units='GHz',/MJy) * kfactor(alpha1,T1) )
   endif else begin
      ; TWO-COMPONENT MODEL: Evaluate equn (6) from Finkbeiner et al

      ; Compute ln(T2) from ln(Rmap).  Loop through the coeffs to save memory.
      lnT2 = fltarr(N_elements(lnR))
      for ii=0, N_elements(RfitA)-1 do lnT2 = lnT2 + RfitA[ii] * lnR^ii
lnR = 0 ; clear memory
      T2 = exp(lnT2)
lnT2 = 0 ; clear memory

      ; Compute T1 as a function of T2; equn (13) of Finkbeiner et al.
      iz1 = (where(Zindx EQ alpha1))[0]
      iz2 = (where(Zindx EQ alpha2))[0]
      h_Pl = 6.6261e-27 ; cm^2 g s^-1
      k_B = 1.3806e-16  ; erg K^-1
      tcoeff = ( (Zintegral[iz2] / (q1q2*Zintegral[iz1])) $
       * (h_Pl*nu100*1.e+9/k_B)^(alpha1-alpha2) )^(1./(4.+alpha1))
      T1 = tcoeff * T2^((4+alpha2)/(4+alpha1))

      Inu = I100 * $
       (f1 * q1q2 * (nu/nu100)^alpha1 * djs_planck(T1,nu,units='GHz',/MJy) $
         + (1-f1) * (nu/nu100)^alpha2 * djs_planck(T2,nu,units='GHz',/MJy) ) / $
       (f1 * q1q2 * djs_planck(T1,nu100,units='GHz',/MJy) * kfactor(alpha1,T1) $
       + (1-f1) * djs_planck(T2,nu100,units='GHz',/MJy) * kfactor(alpha2,T2) )
   endelse
T1 = 0 ; clear memory
T2 = 0 ; clear memory
I100 = 0 ; clear memory

   ; Convert output to requested units
   case strupcase(units) of
      'MJY'   :
      'MICROK': Inu = Inu * fac_flux2temp(nu)
      'THERMO': Inu = Inu * fac_flux2temp(nu) * planckcorr(nu)
      else: begin
          print, "Units must be one of 'MJy', 'microK', or 'thermo'"
          print, '----->  Defaulting to MJy/sr'
      end
   endcase

   ; If OUTFILE is defined, then write to output file
   if (keyword_set(outfile)) then begin
      get_lun, olun
      openw, olun, outfile
      printf, olun, format='(a8,a8,a13,a13)', $
       'l(deg)', 'b(deg)', 'nu(GHz)', 'Inu('+units+')'
      printf, olun, ' ------- ------- ------------ ------------'
      for ii=0, N_elements(gall)-1 do begin
         if (N_elements(nu) EQ 1) then begin
            printf, olun, format='(f8.3,f8.3,e13.5,e13.5)', $
             gall[ii], galb[ii], nu[0], Inu[ii]
         endif else begin
            printf, olun, format='(f8.3,f8.3,e13.5,e13.5)', $
             gall[ii], galb[ii], nu[ii], Inu[ii]
         endelse
      endfor
      close, olun
      free_lun, olun
   endif

   return, Inu
end
;------------------------------------------------------------------------------
