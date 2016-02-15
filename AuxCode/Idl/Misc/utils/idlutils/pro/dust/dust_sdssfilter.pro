;------------------------------------------------------------------------------
;+
; NAME:
;   dust_sdssfilter
;
; PURPOSE:
;   Integrate the extinction curve over a selection of source fuctions and
;   the SDSS filters.
;
; CALLING SEQUENCE:
;   aval = dust_sdssfilter( ixval, [ source=, zsource=, rv=, anorm=, dlam=, /debug ] )
;
; INPUTS:
;   ixval      - Temperature-corrected 100-micron flux (IX) from SFD maps,
;                in MJy/sr [NDATA].
;
; OPTIONAL INPUTS:
;   source     - Source function.  Options are:
;                  'Fstar': Hot F-star model from Kurucz at T_eff=7000K
;                  'Galaxy': SDSS mean galaxy spectrum
;                  'QSO': SDSS mean QSO spectrum
;                Default to 'Galaxy'.
;   zsource    - Redshift of source; default to 0.
;   rv         - Extinction curve parameter R_V; default to 3.1.
;   anorm      - Normalization factor for multiplying IX to obtain
;                extinction at 1 micron; default to SFD normalization
;                of (1.319)*(0.0184) mag/(MJy/sr).
;   dlam       - Spacing of numeric integration in Angstroms; default
;                to 1 Ang.
;   debug      - If set, then make debugging plots with SPLOT.
;   old        - If set, then use the original SDSS filters curves
;                used for the SFD paper.
;
; OUTPUTS:
;   avec       - Extinction in magnitudes [5,NDATA].
;
; COMMENTS:
;
; EXAMPLES:
;   Given a set of Galactic coordinates L,B, evaluate the SFD-predicted
;   extinction in the 5 SDSS filters for a source with the spectrum of an F star:
;     IDL> aval = dust_sdssfilter(dust_getval(l, b, map='IX', /interp), source='Fstar')
;
;   Compare the predicted extinction at a given dust value of IX=10. for two
;   different values of R_V (5.5 vs. 3.1):
;     IDL> print, dust_sdssfilter(10., rv=5.5) / dust_sdssfilter(10.)
;
; PROCEDURES CALLED:
;   dust_intfilter()
;   mrdfits()
;   sxpar()
;
; DATA FILES:
;   $IDLUTILS_DIR/data/filters/kurucz_fstar.dat
;   $IDLUTILS_DIR/data/filters/kpno_atmos.dat
;   $IDLUTILS_DIR/data/filters/sdss_1994_u_noatm.dat
;   $IDLUTILS_DIR/data/filters/sdss_1994_g_noatm.dat
;   $IDLUTILS_DIR/data/filters/sdss_1994_r_noatm.dat
;   $IDLUTILS_DIR/data/filters/sdss_1994_i_noatm.dat
;   $IDLUTILS_DIR/data/filters/sdss_1994_z_noatm.dat
;   $IDLUTILS_DIR/data/filters/sdss_jun2001_u_atm.dat
;   $IDLUTILS_DIR/data/filters/sdss_jun2001_g_atm.dat
;   $IDLUTILS_DIR/data/filters/sdss_jun2001_r_atm.dat
;   $IDLUTILS_DIR/data/filters/sdss_jun2001_i_atm.dat
;   $IDLUTILS_DIR/data/filters/sdss_jun2001_z_atm.dat
;   $IDLUTILS_DIR/data/filters/sdss_meangalaxy_52223.dat
;
; DATA FILES:
;
; REVISION HISTORY:
;   01-Dec-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
function dust_sdssfilter, ixval, source=source, zsource=zsource1, $
 rv=rv, anorm=anorm, dlam=dlam1, debug=debug, old=old

   if (n_params() LT 1) then $
    message, 'Must specify IXVAL'
   if (NOT keyword_set(source)) then source = 'Galaxy'

   ;----------
   ; Generate the list of (5) SDSS filter files

   filtnames = ['u','g','r','i','z']
   if (keyword_set(old)) then begin
      ffiles = filepath('sdss_1994_'+filtnames+'_noatm.dat', $
       root_dir=getenv('IDLUTILS_DIR'), subdirectory=['data','filters'])
      atmfilename1 = filepath('kpno_atmos.dat', $
       root_dir=getenv('IDLUTILS_DIR'), subdirectory=['data','filters'])
   endif else begin
      ffiles = 'sdss_jun2001_'+filtnames+'_atm.dat'
      ffiles = filepath(ffiles, $
       root_dir=getenv('IDLUTILS_DIR'), subdirectory=['data','filters'])
   endelse
   nfile = n_elements(ffiles)

   ;----------
   ; Read the source function

   case source of
   'Fstar': sfilename = filepath('kurucz_fstar.dat', $
              root_dir=getenv('IDLUTILS_DIR'), subdirectory=['data','filters'])
   'Galaxy': sfilename = filepath('sdss_meangalaxy_52223.dat', $
              root_dir=getenv('IDLUTILS_DIR'), subdirectory=['data','filters'])
   'QSO': sfilename = filepath('sdss_meanqso_52223.dat', $
              root_dir=getenv('IDLUTILS_DIR'), subdirectory=['data','filters'])
   else: message, 'Unknown SOURCE'
   endcase

   ;----------
   ; Loop over each filter and compute the extinction

   ndat = n_elements(ixval)
   aval = fltarr(nfile, ndat)
   for ifile=0, nfile-1 do begin
      f_tabulate = 0
      s_tabulate = 0
      aval[ifile,*] = dust_intfilter(ixval, ffilename=ffiles[ifile], $
       sfilename=sfilename, atmfilename=atmfilename1, $
       rv=rv, anorm=anorm, dlam=dlam1, zsource=zsource1, debug=debug)
   endfor

   return, aval
end
;------------------------------------------------------------------------------
