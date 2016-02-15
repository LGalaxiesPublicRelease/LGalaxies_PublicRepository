;------------------------------------------------------------------------------
;+
; NAME:
;   tycho_read
;
; PURPOSE:
;   Read the Tycho catalog (and generate IDL save sets)
;
; CALLING SEQUENCE:
;   bigdat = tycho_read( [ racen=, deccen=, radius=, node=, incl=, hwidth=, $
;    epoch=, columns= ] )
;
; OPTIONAL INPUTS:
;   racen:       Central RA for selecting a region of stars [J2000 deg]
;   deccen:      Central DEC for selecting a region of stars [J2000 deg]
;   radius:      Radius for selecting a region of stars [deg]
;   node:        Node of great circle for selecting a stripe of stars [deg]
;   incl:        Inclination of great circle for selecting a stripe [deg]
;   hwidth:      Half-width of great circle for selecting a stripe [deg]
;   epoch:       If set, then apply proper motion correction from epoch 2000.
;                to that set by EPOCH.  Note that equinox is always J2000.
;   columns:     Select only these columns from the data files.
;                Must include 'ramdeg' and 'demdeg' if a circular or
;                great circle region is selected.
;                Must also include 'pmra' and 'pmde' if EPOCH is specified.
;
; OUTPUTS:
;   bigdat:      Returned structure containing data
;
; COMMENTS:
;   To trim to objects within a circle, RACEN, DECCEN and RADIUS must
;   all be set.
;
;   To trim to objects along a great circle, NODE, INCL and HWIDTH must
;   all be set.
;
;   The equinox of the returned data is always J2000, even if the epoch
;   has been chosen to be different.
;
;   Note that 104,189 of the stars (4.3%) do not have a mean position,
;   and RAMDEG,DEMDEG are set to exactly zero.  These stars still have
;   an observed position in the RADEG,DEDEG fields, which are typically
;   from epoch ~1991.  The expected proper motion from that epoch to
;   2000 is approximately 0.08 arcsec for stars in this catalog, or
;   about 0.01 arcsec/yr.
;
; BUGS:
;
; DATA FILES:
;   The following data files can be copied from:
;     http://adc.gsfc.nasa.gov/adc-cgi/cat.pl?/catalogs/1/1259/
;   and should be put in a directory pointed to by the environment
;   variable $TYCHO2_DIR:
;     $TYCHO2_DIR/index.dat.gz
;     $TYCHO2_DIR/tyc2_*.dat.gz
;   These files may be either left compressed, or uncompressed.
;   This code looks for either.
;
; INTERNAL SUPPORT ROUTINES:
;   tyc_convert()
;   tyc_rdindex
;   tyc_append_list
;   tyc_readascii()
;   tyc_readfits()
;
; PROCEDURES CALLED:
;   djs_diff_angle()
;   djs_filepath()
;   hip_epoch
;   mrdfits()
;   mwrfits
;   radec_to_munu
;   splog
;   tycho_epoch
;
; REVISION HISTORY:
;   Written D. Schlegel, 31 December 2002, Princeton
;-
;------------------------------------------------------------------------------
function tyc_convert, barr, format

   stringval = string(barr)
   ibad = where(strcompress(/remove_all,stringval) EQ '')

   case format of
   'I': begin
      if (ibad[0] NE -1) then stringval[ibad] = '0'
      retval = long(stringval)
      end
   'F': begin
      if (ibad[0] NE -1) then stringval[ibad] = '0'
      retval = float(stringval)
      end
   'D': begin
      if (ibad[0] NE -1) then stringval[ibad] = '0'
      retval = double(stringval)
      end
   'A': begin
      retval = string(barr)
      end
   endcase

   return, retval
end
;------------------------------------------------------------------------------
pro tyc_rdindex, tyc_dir=tyc_dir

   common com_tyc_read, tycindex, ra_block, dec_block, rad_block

   ;----------
   ; Find the Tycho index file.  First look for the g-zipped version,
   ; if not found then look for the uncompressed version

   indfile = (findfile(filepath('index.dat.gz', root_dir=tyc_dir)))[0]
   compress = 1
   if (NOT keyword_set(indfile)) then begin
      indfile = (findfile(filepath('index.dat', root_dir=tyc_dir)))[0]
      compress = 0
   endif
   if (NOT keyword_set(indfile)) then begin
      splog, 'WARNING: Tycho-2 index file not found.'
      splog, 'You may download these data from:'
      splog, '  http://adc.gsfc.nasa.gov/adc-cgi/cat.pl?/catalogs/1/1259/'
      return
   endif

   ;----------
   ; Read the entire Tycho index file into a byte array

;   nline = numlines(indfile) ; NUMLINES fails if the index file is compressed
   nline = 9538
;   splog, 'Reading Tycho index file ' + indfile
   barr = bytarr(43,nline)
   tmpString = ''
   get_lun, ilun
   openr, ilun, indfile, compress=compress
   for iline=0L, nline-1 do begin
      readf, ilun, tmpString
      barr[0:strlen(tmpString)-1,iline] = byte( tmpString )
   endfor
   close, ilun
   free_lun, ilun

   ; Create data structure and allocate memory for data.
   pTemp = $
   { tycind, $
     rec_t2    : 0L            , $  ; Tycho-2 record of 1st star in region
     rec_s1    : 0L            , $  ; Suppl-1 record of 1st star in region
     RAmin     : 0.0           , $  ; Smallest RA in region
     RAmax     : 0.0           , $  ; Largest RA in region
     DEmin     : 0.0           , $  ; Smallest DEC in region
     DEmax     : 0.0             $  ; Largest DEC in region
   }
   tycindex = replicate( {tycind}, nline )

   ; Copy byte data into data structure
   tycindex.rec_t2   = tyc_convert(barr[  0:  6,*], 'I')
   tycindex.rec_s1   = tyc_convert(barr[  8: 13,*], 'I')
   tycindex.RAmin    = tyc_convert(barr[ 15: 20,*], 'F')
   tycindex.RAmax    = tyc_convert(barr[ 22: 27,*], 'F')
   tycindex.DEmin    = tyc_convert(barr[ 29: 34,*], 'F')
   tycindex.DEmax    = tyc_convert(barr[ 36: 41,*], 'F')

   ra_block = 0.5 * (tycindex.ramin + tycindex.ramax)
   dec_block = 0.5 * (tycindex.demin + tycindex.demax)
   ang1 = djs_diff_angle(tycindex.ramin, tycindex.demin, $
    tycindex.ramax, tycindex.demax)
   ang2 = djs_diff_angle(tycindex.ramin, tycindex.demax, $
    tycindex.ramax, tycindex.demin)
   rad_block = 0.5 * (ang1 > ang2)

   return
end
;------------------------------------------------------------------------------
pro tyc_append_list, filelist, filenum, row1, row2

   if (keyword_set(filelist)) then nlist = n_elements(filelist) $
    else nlist = 0
   if (nlist GT 0) then $
    filelist = [filelist, filelist[0]] $
   else $
    filelist = create_struct('filename', '', 'row1', 0L, 'row2', 0L)
   filelist[nlist].filename = string(filenum, format='("tyc2_",i2.2)')
   filelist[nlist].row1 = row1
   filelist[nlist].row2 = row2

   return
end
;------------------------------------------------------------------------------
function tyc_readascii, filename, range=range1, columns=columns

   ;----------
   ; Now look for either a compressed or un-compressed ASCII version.

   thisfile = (findfile(filename+'.dat.gz'))[0]
   compress = 1
   if (NOT keyword_set(thisfile)) then begin
      thisfile = (findfile(filename+'.dat'))[0]
      compress = 0
   endif
   if (NOT keyword_set(thisfile)) then begin
      splog, 'WARNING: Tycho-2 data file not found ' + filename+'.dat .'
      splog, 'You may download these data from:'
      splog, '  http://adc.gsfc.nasa.gov/adc-cgi/cat.pl?/catalogs/1/1259/'
      return, 0
   endif

   ;----------
   ; Read a single Tycho data file into a byte array

   if (thisfile NE 'tyc2_19') then expectlines = 127000L $
    else expectlines = 126913L
;  NUMLINES fails if the index file is compressed
;   if (keyword_set(range1)) then range = range1 $
;    else range = [0L, numlines(thisfile)-1]
   if (keyword_set(range1)) then range = range1 $
    else range = [0L, expectlines-1]
   nline = range[1] - range[0] + 1
;   splog, 'Reading Tycho data file ' + thisfile + ' RANGE=', $
;    range[0], range[1]
   barr = bytarr(207,nline)
   tmpString = ''
   get_lun, ilun
   openr, ilun, thisfile
   point_lun, ilun, 207L * range[0]
   for iline=0L, nline-1 do begin
      readf, ilun, tmpString
      barr[0:strlen(tmpString)-1,iline] = byte( tmpString )
   endfor
   close, ilun
   free_lun, ilun

   ; Create data structure and allocate memory for data.
   pTemp = $
   { tycobj, $
     TYC1      : 0L            , $  ; TYC1 from TYC or GSC [1,9537]
     TYC2      : 0L            , $  ; TYC2 from TYC or GSC [1,12121]
     TYC3      : 0             , $  ; TYC3 from TYC [1,3]
     pflag     : ' '           , $  ; Mean position flag [ PX]
     RAmdeg    : 0.d0          , $  ; Mean RA, ICRS, epoch=J2000 [deg]
     DEmdeg    : 0.d0          , $  ; Mean DEC, ICRS, epoch=J2000 [deg]
     pmRA      : 0.0           , $  ; Prop. mot. in RA*cos(dec) [mas/yr]
     pmDE      : 0.0           , $  ; Prop. mot. in DEC [mas/yr]
     e_RAmdeg  : 0             , $  ; Err of RA*cos(dec) at mean epoch [mas]
     e_DEmdeg  : 0             , $  ; Err of DEC at mean epoch [mas]
     e_pmRA    : 0.0           , $  ; Err of prop. mot. in RA*cos(dec) [mas/yr]
     e_pmDE    : 0.0           , $  ; Err of prop. mot. in DEC [mas/yr]
     EpRAm     : 0.0           , $  ; Mean epoch of RA
     EpDEm     : 0.0           , $  ; Mean epoch of DEC
     Num       : 0             , $  ; Number of positions used [2,36]
     q_RAmdeg  : 0.0           , $  ; Goodness of fit for mean RA [0,9.9]
     q_DEmdeg  : 0.0           , $  ; Goodness of fit for mean DEC [0,9.9]
     q_pmRA    : 0.0           , $  ; Goodness of fit for pmRA [0,9.9]
     q_pmDE    : 0.0           , $  ; Goodness of fit for pmDE [0,9.9]
     BTmag     : 0.0           , $  ; Tycho-2 BT magnitude
     e_BTmag   : 0.0           , $  ; Err of BTmag
     VTmag     : 0.0           , $  ; Tycho-2 VT magnitude
     e_VTmag   : 0.0           , $  ; Err of VTmag
     prox      : 0             , $  ; Proximity indicator [3,999]
     TYC       : ' '           , $  ; Tycho-1 star [ T]
     HIP       : 0L            , $  ; Hipparcos number [1,120404]
     CCDM      : '   '         , $  ; CCDM component identifier for HIP stars
     RAdeg     : 0.d0          , $  ; Observed Tycho-2 RA, ICRS [deg]
     DEdeg     : 0.d0          , $  ; Observed Tycho-2 DEC, ICRS [deg]
     EpRA_1990 : 0.0           , $  ; Epoch-1990 of RAdeg [0.81,2.13]
     EpDE_1990 : 0.0           , $  ; Epoch-1990 of DEdeg [0.72,2.36]
     e_RAdeg   : 0.0           , $  ; Err of observed Tycho-2 RA*cos(dec) [mas]
     e_DEdeg   : 0.0           , $  ; Err of observed Tycho-2 DEC [mas]
     posflg    : ' '           , $  ; Type of Tycho-2 solution [ DP]
     corr      : ' '             $  ; Correlation (RAdeg,DEdeg) [-1,1]
   }
   bigdat = replicate( {tycobj}, nline )

   ; Copy byte data into data structure
   bigdat.TYC1     = tyc_convert(barr[  0:  3,*], 'I')
   bigdat.TYC2     = tyc_convert(barr[  5:  9,*], 'I')
   bigdat.TYC3     = tyc_convert(barr[ 11: 11,*], 'I')
   bigdat.pflag    = tyc_convert(barr[ 13: 13,*], 'A')
   bigdat.RAmdeg   = tyc_convert(barr[ 15: 26,*], 'D')
   bigdat.DEmdeg   = tyc_convert(barr[ 28: 39,*], 'D')
   bigdat.pmRA     = tyc_convert(barr[ 41: 47,*], 'F')
   bigdat.pmDE     = tyc_convert(barr[ 49: 55,*], 'F')
   bigdat.e_RAmdeg = tyc_convert(barr[ 57: 59,*], 'I')
   bigdat.e_DEmdeg = tyc_convert(barr[ 61: 63,*], 'I')
   bigdat.e_pmRA   = tyc_convert(barr[ 65: 68,*], 'F')
   bigdat.e_pmDE   = tyc_convert(barr[ 70: 73,*], 'F')
   bigdat.EpRAm    = tyc_convert(barr[ 75: 81,*], 'F')
   bigdat.EpDEm    = tyc_convert(barr[ 83: 89,*], 'F')
   bigdat.Num      = tyc_convert(barr[ 91: 92,*], 'I')
   bigdat.q_RAmdeg = tyc_convert(barr[ 94: 96,*], 'F')
   bigdat.q_DEmdeg = tyc_convert(barr[ 98:100,*], 'F')
   bigdat.q_pmRA   = tyc_convert(barr[102:104,*], 'F')
   bigdat.q_pmDE   = tyc_convert(barr[106:108,*], 'F')
   bigdat.BTmag    = tyc_convert(barr[110:115,*], 'F')
   bigdat.e_BTmag  = tyc_convert(barr[117:121,*], 'F')
   bigdat.VTmag    = tyc_convert(barr[123:128,*], 'F')
   bigdat.e_VTmag  = tyc_convert(barr[130:134,*], 'F')
   bigdat.prox     = tyc_convert(barr[136:138,*], 'I')
   bigdat.TYC      = tyc_convert(barr[140    ,*], 'A')
   bigdat.HIP      = tyc_convert(barr[142:147,*], 'I')
   bigdat.CCDM     = tyc_convert(barr[148:150,*], 'A')
   bigdat.RAdeg    = tyc_convert(barr[152:163,*], 'D')
   bigdat.DEdeg    = tyc_convert(barr[165:176,*], 'D')
   bigdat.EpRA_1990= tyc_convert(barr[178:181,*], 'F')
   bigdat.EpDE_1990= tyc_convert(barr[183:186,*], 'F')
   bigdat.e_RAdeg  = tyc_convert(barr[188:192,*], 'F')
   bigdat.e_DEdeg  = tyc_convert(barr[194:198,*], 'F')
   bigdat.posflg   = tyc_convert(barr[200    ,*], 'A')
   bigdat.corr     = tyc_convert(barr[202:205,*], 'F')

   ; Delete the byte data array
   barr = 0

   ;----------
   ; If COLUMNS is set, then trim the structure to the specified columns

   if (keyword_set(columns)) then begin
      tags = tag_names(bigdat)
      ntag = n_elements(tags)
      qtag = bytarr(ntag)
      for icol=0, n_elements(columns)-1 do $
       qtag = qtag OR tags EQ strupcase(columns[icol])
      ikeep = where(qtag, nkeep)
      if (nkeep EQ 0) then ikeep = lindgen(ntag) ; Keep all columns if no match
      trimdat = create_struct(tags[ikeep[0]], bigdat[0].(ikeep[0]))
      for i=1, nkeep-1 do $
       trimdat = create_struct(trimdat, tags[ikeep[i]], bigdat[0].(ikeep[i]))
      trimdat = replicate(trimdat, n_elements(bigdat))
      struct_assign, bigdat, trimdat
      bigdat = temporary(trimdat)
   endif

   return, bigdat
end
;------------------------------------------------------------------------------
; If a FITS data file exists, return the name of that file.
; Otherwise, generate the FITS file if it is possible to write to disk,
; and return the name of that newly-written FITS file.
; Otherwise, return null string.

function tyc_readfits, filename

   ;----------
   ; Return the name of the FITS file (or compressed file) if it exists.

   thisfile = (findfile(filename+'.fits*'))[0]
   if (keyword_set(thisfile)) then return, thisfile

   ;----------
   ; Test whether we can write to disk.  If not, then return a null string.

   fitsfile = filename+'.fits'
   openw, olun, fitsfile, /get_lun, /delete, error=err
   if (keyword_set(olun)) then close, olun
   if (keyword_set(err)) then begin
      splog, 'Unable to write FITS file to disk for ' + filename
      return, ''
   endif

   ;----------
   ; Read the full ASCII file, and write a FITS file.
   ; Return the name of that FITS file.

   splog, 'Generating FITS version of ' + filename
   thisdat = tyc_readascii(filename)
   if (keyword_set(thisdat)) then begin
      mwrfits, thisdat, fitsfile, /create
   endif else begin
      splog, 'Failure to generate FITS file of ' + filename
      fitsfile = ''
   endelse

   return, fitsfile
end
;------------------------------------------------------------------------------
function tyc_readfile, filename, range=range, columns=columns

   ;----------
   ; Find a single Tycho data file.  First look for a FITS version.
   ; If it does not exist (and cannot be created), then read the ASCII version.

   thisfile = tyc_readfits(filename)
   if (keyword_set(thisfile)) then begin
;      splog, 'Reading FITS file ' + thisfile + ' RANGE=', range[0], range[1]
      bigdat = mrdfits(thisfile, 1, range=range, columns=columns, /silent)
   endif else begin
      bigdat = tyc_readascii(filename, range=range, columns=columns)
   endelse

   return, bigdat
end
;------------------------------------------------------------------------------
function tycho_read, racen=racen, deccen=deccen, radius=radius, $
 node=node, incl=incl, hwidth=hwidth, epoch=epoch, columns=columns

   common com_tyc_read, tycindex, ra_block, dec_block, rad_block

   tyc_dir = getenv('TYCHO2_DIR')
   if (NOT keyword_set(tyc_dir)) then begin
      splog, 'WARNING: TYCHO2_DIR environment variable not set.'
      splog, 'This must point to a data directory with files downloade from:'
      splog, '  http://adc.gsfc.nasa.gov/adc-cgi/cat.pl?/catalogs/1/1259/'
      return, 0
   endif

   ;----------
   ; Read the Tycho index file if it is not already cached in the common block

   if (NOT keyword_set(tycindex)) then $
    tyc_rdindex, tyc_dir=tyc_dir
   if (NOT keyword_set(tycindex)) then return, 0
   nindex = n_elements(tycindex)

   ;----------
   ; Pad the search distances by some amount.
   ; There is 0.01 degrees, since one of the blocks of Tycho data
   ; that should start at RA=0 actually starts at RA=-0.01 deg.
   ; Also pad by the maximum proper motion that could occur,
   ; which is 10.3083 arcsec/yr = 0.002863 deg/yr.

   dpad = 0.01d0
   if (keyword_set(epoch)) then dpad = dpad + 0.00287 * abs(epoch - 2000.)

   ;----------
   ; Determine which blocks of which data files to read

   qkeep = bytarr(nindex) + 1B

   qreadcircle = (n_elements(racen) GT 0 AND n_elements(deccen) GT 0 $
    AND n_elements(radius) GT 0)
   qreadstripe = (n_elements(node) GT 0 AND n_elements(incl) GT 0 $
    AND n_elements(hwidth) GT 0)

   if (qreadcircle) then begin
      decmin = deccen - (radius+dpad)
      decmax = deccen + (radius+dpad)
      if (decmax GE 90) then cosdec = 1.d0 $
       else cosdec = cos(decmax/!radeg)
      ramin = racen - (radius+dpad) / cosdec
      ramax = racen + (radius+dpad) / cosdec
      if (ramax - ramin GE 360 OR decmax GE 90) then begin
         ramin = 0.0
         ramax = 360.0
      endif

      qkeepdec = 1B - (decmin GT tycindex.demax OR decmax LT tycindex.demin)
      qkeepra = 1B - (ramin GT tycindex.ramax OR ramax LT tycindex.ramin)
      if (ramin LT 0) then $
       qkeepra = qkeep $
       OR (1B - (ramin+360 GT tycindex.ramax OR ramax+360 LT tycindex.ramin))
      if (ramax GT 360) then $
       qkeepra = qkeep $
       OR (1B - (ramin-360 GT tycindex.ramax OR ramax-360 LT tycindex.ramin))
      qkeep = qkeep AND qkeepdec AND qkeepra
   endif
   if (qreadstripe) then begin
      radec_to_munu, ra_block, dec_block, node=node, incl=incl, $
       mu, nu
      qkeep = qkeep AND (abs(nu) LE (hwidth + rad_block + dpad))
   endif

   qkeep[nindex-1] = 0 ; The last index block contains no objects
   if (total(qkeep) EQ 0) then return, 0

   ;----------
   ; Convert these blocks of data into contiguous record numbers

   qpadkeep = [0B, qkeep, 0B]
   i1 = where(qpadkeep[1:nindex] GT qpadkeep[0:nindex-1], nblock)
   i2 = where(qpadkeep[1:nindex] GT qpadkeep[2:nindex+1])
   record1 = tycindex[i1].rec_t2 - 1L ; Change from 1-indexed to 0-indexed
   record2 = tycindex[i2+1].rec_t2 - 2L ; Change from 1-indexed to 0-indexed

   ;----------
   ; Convert these record numbers into file names and row numbers

   filelist = 0
   for iblock=0, nblock-1 do begin
      file1 = long(record1[iblock] / 127000.)
      row1 = record1[iblock] MOD 127000L
      file2 = long(record2[iblock] / 127000.)
      row2 = record2[iblock] MOD 127000L
      if (file1 EQ file2) then begin
         tyc_append_list, filelist, file1, row1, row2
      endif else begin
         tyc_append_list, filelist, file1, row1, 126999L
         for ifile=file1+1, file2-1  do $
          tyc_append_list, filelist, ifile, 0L, 126999L
         tyc_append_list, filelist, file2, 0L, row2
      endelse
   endfor

   ;----------
   ; Read the files

   ntot = long(total(filelist.row2 - filelist.row1 + 1L))
   itot = 0L
   nlist = n_elements(filelist)
   for ifile=0, nlist-1 do begin
      thisfile = djs_filepath(filelist[ifile].filename, root_dir=tyc_dir)
      row1 = filelist[ifile].row1
      row2 = filelist[ifile].row2
      thisdat = tyc_readfile(thisfile, range=[row1,row2], columns=columns)
      nthis = n_elements(thisdat)
      if (NOT keyword_set(bigdat)) then begin
         blankdat = thisdat[0]
         struct_assign, {junk:0}, blankdat
         bigdat = replicate(blankdat, ntot)
      endif
      bigdat[itot:itot+nthis-1] = thisdat
      itot = itot + nthis
   endfor
   if (itot NE ntot) then $
    message, 'Error encountered reading Tycho data files'

   ;----------
   ; Transform to a different epoch

   if (keyword_set(epoch)) then tycho_epoch, epoch, bigdat

   ;----------
   ; Further trim to the requested circle on the sky

   if (qreadcircle OR qreadstripe) then begin
      qkeep = bytarr(ntot) + 1B
      if (qreadcircle) then begin
         adiff = djs_diff_angle(bigdat.ramdeg, bigdat.demdeg, $
          racen, deccen)
         qkeep = qkeep AND (adiff LE radius)
      endif
      if (qreadstripe) then begin
         radec_to_munu, bigdat.ramdeg, bigdat.demdeg, node=node, incl=incl, $
          mu, nu
         qkeep = qkeep AND (abs(nu) LE hwidth)
      endif
      ikeep = where(qkeep, nkeep)
      if (nkeep EQ 0) then bigdat = 0 $
       else bigdat = bigdat[ikeep]
   endif

   return, bigdat
end
;------------------------------------------------------------------------------
