;------------------------------------------------------------------------------
;+
; NAME:
;   predict_sync
;
; PURPOSE:
;   Read predictions of synchtron emission and return sky intensity.
;
; CALLING SEQUENCE:
;   value = predict_sync( [ gall, galb, nu=nu, infile=infile, $
;    skipline=skipline, outfile=outfile, interp=interp, $
;    noloop=noloop, verbose=verbose, ipath=ipath, units=units ] )
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   gall:       Galactic longitude(s) in degrees
;   galb:       Galactic latitude(s) in degrees
;   nu:         Frequency in GHz.  If this is a vector, it must be the same
;               dimension as GALL and GALB.  If this is a scalar, then it
;               applies to all positions GALL, GALB.
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
;   ipath:      Path name for synchrotron maps; default to path set by the
;               environment variable $DUST_DIR/map, or to the current
;               directory.
;   units:      Units for output values:
;               'MJy'    : MJy/sr (default)
;               'microK' : brightness (antenna) temperature [micro-Kelvin]
;               'thermo' : thermodynamic temperature [micro-Kelvin]
;
; OUTPUTS:
;   value:      Predicted emission value(s) from synchtrotron maps
;
; COMMENTS:
;   These predictions are based upon the following radio surveys:
;     408 MHz - Haslam et al. 1981, A&A, 100, 209
;     1420 MHz - Reich & Reich 1986, A&AS, 63, 205
;     2326 MHz - Jonas, Baart, & Nicolson 1998 MNRAS, 297, 977
;   A detailed description of this data product may be found in:
;   "Extrapolation of the Haslam 408 MHz survey to CMBR Frequencies",
;   by D. P. Finkbeiner, & M. Davis, in preparation.
;
;   Either the coordinates GALL, GALB and NU must be set, or their values
;   must exist in the file INFILE.  Output is written to the variable VALUE
;   and/or the file OUTFILE.
;
; EXAMPLES:
;   Read the predicted synchtrotron emission at Galactic (l,b)=(12,+34.5) 
;   interpolating from the nearest 4 pixels, and return result in VALUE.
;   > value = predict_sync(12, 34.5, nu=3000, /interp)
;
; PROCEDURES CALLED:
;   fac_flux2temp()
;   planckcorr()
;   rdfloat
;   wcs_getval()
;
; DATA FILES:
;   Haslam_clean_ngp.fits
;   Haslam_clean_sgp.fits
;   Synch_Beta_ngp.fits
;   Synch_Beta_sgp.fits
;
; REVISION HISTORY:
;   17-Mar-1999  Written by David Schlegel, Princeton
;                & Doug Finkbeiner, Berkeley
;-
;------------------------------------------------------------------------------
function predict_sync, gall, galb, nu=nu, infile=infile, $
 skipline=skipline, outfile=outfile, resolution=resolution, model=model, $
 interp=interp, noloop=noloop, verbose=verbose, ipath=ipath, units=units

   if (NOT keyword_set(infile) $
    AND (NOT keyword_set(gall) OR NOT keyword_set(galb) $
         OR NOT keyword_set(nu)) $
    ) then begin
      print, 'Must set either coordinates and frequency or INFILE'
      return, -1
   endif
   if (NOT keyword_set(units)) then units = 'MJy'

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

   if (NOT keyword_set(ipath)) then begin      dust_dir = getenv('DUST_DIR')
      if (dust_dir NE '') then ipath = dust_dir+'/maps/' $
       else ipath = './'
   endif

   ; Read the Haslam map (units of K brightness temp) this Haslam has had
   ; 2.73K subtracted for the CMB (no dipole), point sources removed, and
   ; it is Fourier destriped. 

   Amap = wcs_getval( $
    ['Haslam_clean_ngp.fits', 'Haslam_clean_sgp.fits'], $
    gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)
   Bmap = wcs_getval( $
    ['Synch_Beta_ngp.fits', 'Synch_Beta_sgp.fits'], $
    gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)

   ; microK brightness temp (Beta map is actually negative of spectral
   ; index, and ranges from roughly 2.5 < beta < 3.0)

   Inu = 1E6 * Amap * (0.408/nu)^Bmap

   ; Convert output to requested units
   case units of
      'MJy'   : Inu = Inu / fac_flux2temp(nu)
      'microK': Inu = Inu
      'thermo': Inu = Inu * planckcorr(nu)
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
