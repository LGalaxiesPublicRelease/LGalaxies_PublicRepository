;+
; NAME:
;   ucac_readzone()
;
; PURPOSE:
;   Read the raw UCAC data files for a specific declination zone within
;   a given RA range.
;
; CALLING SEQUENCE:
;   outdat = ucac_readzone(zone, ra_min, ra_max)
;
; INPUTS:
;   zone       - UCAC zone number (corresponding to a particular declination)
;   ra_min     - Minimum RA [deg]
;   ra_max     - Maximum RA [deg]
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;   outdat     - Structure with UCAC data in its raw catalog format;
;                return 0 if no stars found
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Quantities are de-coded to meaningful units, e.g. converting RA to degrees.
;
; EXAMPLES:
;
; BUGS:
;   The values of TWOMASS_PH,TWOMASS_CC do not appear to make any
;   sense relative to the UCAC documentation.
;
; PROCEDURES CALLED:
;   ucac_readindex()
;
; REVISION HISTORY:
;   27-May-2003  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function ucac_readzone, thiszone, ra_min, ra_max

   common com_ucac, uindex

   ;----------
   ; Check inputs

   if (n_params() LT 3) then begin
      print, 'Wrong number of parameters!'
      return, 0
   endif
   if (ra_min GT ra_max OR ra_min LT 0 OR ra_min GT 360 $
    OR ra_max LT 0 OR ra_max GT 360) then begin
      print, 'Invalid RA_MIN,RA_MAX'
      return, 0
   endif

   ucac_dir = getenv('UCAC_DIR')
   if (NOT keyword_set(ucac_dir)) then begin
      print, 'Environment variable UCAC_DIR must be set!'
      return, 0
   endif

   ;----------
   ; Read the index file

   uindex = ucac_readindex()

   ;----------
   ; Determine where to seek in this zone file.

   jj = where(uindex.zn EQ thiszone, ct)
   if (ct EQ 0) then begin
      print, 'This zone not found'
      return, 0
   endif

   j1 = (where(uindex[jj].ramax * 15.d GE ra_min))[0]
   j1 = j1 > 0L
   j2 = (reverse(where(uindex[jj].ramax * 15.d LE ra_max)))[0] + 1L
   j2 = j2 < (n_elements(jj) - 1L) ; In the case that RA_MAX=360 deg

   if (j1 EQ 0) then i1 = 0L $
    else i1 = uindex[jj[j1-1]].naz
   i2 = uindex[jj[j2]].naz - 1L
   nrecord = i2 - i1 + 1L
   if (nrecord EQ 0) then return, 0

   ;----------
   ; Read the binary format data

   thisfile = filepath(string(thiszone,format='("z",i3.3)'), $
    root_dir=ucac_dir, subdir='u2')

   blankdat = create_struct( $
    'RA'    , 0L, $
    'DEC'   , 0L, $
    'RMAG'  , 0 , $
    'E_RAM' , 0B, $
    'E_DEM' , 0B, $
    'NOBS'  , 0B, $
    'RFLAG' , 0B, $
    'NCAT'  , 0B, $
    'CFLAG' , 0B, $
    'EPRAM' , 0 , $
    'EPDEM' , 0 , $
    'PMRA'  , 0L, $
    'PMDE'  , 0L, $
    'E_PMRA', 0B, $
    'E_PMDE', 0B, $
    'Q_PMRA', 0B, $
    'Q_PMDE', 0B, $
    'TWOMASS_ID' , 0L, $
    'TWOMASS_J'  , 0 , $
    'TWOMASS_H'  , 0 , $
    'TWOMASS_KS' , 0 , $
    'TWOMASS_PH' , 0B, $
    'TWOMASS_CC' , 0B)
   rawdat = replicate(blankdat, nrecord)
   openr, ilun, thisfile, /get_lun, /swap_if_big_endian
   point_lun, ilun, i1 * n_tags(blankdat, /length)
   readu, ilun, rawdat
   close, ilun
   free_lun, ilun

   ;----------
   ; Convert to rational units

   outdat = replicate( create_struct( $
    'RAMDEG' , 0d, $
    'DEMDEG' , 0d, $
    'RMAG'   , 0., $
    'E_RAM'  , 0., $
    'E_DEM'  , 0., $
    'NOBS'   , 0 , $
    'RFLAG'  , 0B, $
    'NCAT'   , 0 , $
    'CFLAG'  , 0B, $
    'EPRAM'  , 0., $
    'EPDEM'  , 0., $
    'PMRA'   , 0., $
    'PMDE'   , 0., $
    'E_PMRA' , 0., $
    'E_PMDE' , 0., $
    'Q_PMRA' , 0., $
    'Q_PMDE' , 0., $
    'TWOMASS_ID' , 0L, $
    'TWOMASS_J'  , 0., $
    'TWOMASS_H'  , 0., $
    'TWOMASS_KS' , 0., $
    'TWOMASS_PH' , 0, $
    'TWOMASS_CC' , 0), nrecord)

   outdat.ramdeg = rawdat.ra / (3600d0*1000d0) ; Mean RA, ICRS, epoch=J2000 [deg]
   outdat.demdeg = rawdat.dec / (3600d0*1000d0) ; Mean DEC, ICRS, epoch=J2000 [deg]
   outdat.rmag = rawdat.rmag / 100. ; R-band magnitude
   outdat.e_ram = rawdat.e_ram - 127. ; Err of RA*cos(dec) at mean epoch [mas]
   outdat.e_dem = rawdat.e_dem - 127. ; Err of DEC at mean epoch [mas]
   outdat.nobs = rawdat.nobs
   outdat.rflag = rawdat.rflag
   outdat.ncat = rawdat.ncat
   outdat.cflag = rawdat.cflag
   outdat.epram = rawdat.epram / 1000. + 1975. ; Mean epoch of RA [yr]
   outdat.epdem = rawdat.epdem / 1000. + 1975. ; Mean epoch of DEC [yr]
   outdat.pmra = rawdat.pmra / 10. ; Proper motion in RA (no cos(dec)) [mas/yr]
   outdat.pmde = rawdat.pmde / 10. ; Proper motion in DEC [mas/yr]
   outdat.e_pmra = (rawdat.e_pmra - 127.) / 10. ; Error of prop. mot. in RA*cosd [mas/yr]
   outdat.e_pmde = (rawdat.e_pmde - 127.) / 10. ; Error of prop. mot. in DEC [mas/yr]
   outdat.q_pmra = (rawdat.q_pmra - 127.) / 20. ; Goodness of fit for pmra
   outdat.q_pmde = (rawdat.q_pmde - 127.) / 20. ; Goodness of fit for pmde
   outdat.twomass_id = rawdat.twomass_id ; 2MASS pts_key star identifier
   outdat.twomass_j = rawdat.twomass_j / 1000. ; 2MASS J magnitude
   outdat.twomass_h = rawdat.twomass_h / 1000. ; 2MASS H magnitude
   outdat.twomass_ks = rawdat.twomass_ks / 1000. ; 2MASS Ks magnitude
   outdat.twomass_ph = rawdat.twomass_ph - 127 ; 2MASS modified ph_qual flag
   outdat.twomass_cc = rawdat.twomass_cc - 127 ; 2MASS modified cc_flg

   ;----------
   ; Trim to the RA range requested

   ikeep = where(outdat.ramdeg GE ra_min AND outdat.ramdeg LE ra_max, nkeep)
   if (nkeep EQ 0) then begin
      return, 0
   endif
   outdat = outdat[ikeep]

   return, outdat
end
;------------------------------------------------------------------------------
