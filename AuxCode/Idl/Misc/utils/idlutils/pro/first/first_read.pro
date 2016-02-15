;------------------------------------------------------------------------------
;+
; NAME:
;   first_read
;
; PURPOSE:
;   Read the FIRST catalog (and generate IDL save sets)
;
; CALLING SEQUENCE:
;   bigdat = first_read( [ racen, deccen, radius, node=, incl=, hwidth=, $
;    columns= ] )
;
; OPTIONAL INPUTS:
;   racen:       Central RA for selecting a region of stars [J2000 deg]
;   deccen:      Central DEC for selecting a region of stars [J2000 deg]
;   radius:      Radius for selecting a region of stars [deg]
;   node:        Node of great circle for selecting a stripe of stars [deg]
;   incl:        Inclination of great circle for selecting a stripe [deg]
;   hwidth:      Half-width of great circle for selecting a stripe [deg]
;   columns:     Select only these columns from the data files.
;                Must include 'FIRST_RA' and 'FIRST_DEC' if a circular or
;                great circle region is selected.
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
;   The equinox of the returned data is always J2000.
;
; BUGS:
;
; DATA FILES:
;   The following data files can be copied from:
;     http://sundog.stsci.edu/first/catalogs/
;   and should be put in a directory pointed to by the environment
;   variable $FIRST_DIR:
;     $FIRST_DIR/catalog_03apr11.bin
;   This file must be uncompressed.
;
; INTERNAL SUPPORT ROUTINES:
;   first_convert()
;   first_testwrite()
;   first_readascii()
;   first_rdindex
;   first_readfits()
;   first_readfile()
;
; PROCEDURES CALLED:
;   djs_diff_angle()
;   mrdfits()
;   mwrfits
;   radec_to_munu
;   splog
;
; REVISION HISTORY:
;   Written D. Schlegel, 18 July 2003, Princeton
;-
;------------------------------------------------------------------------------
function first_convert, barr, format

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
function first_testwrite, filename

   ;----------
   ; Test whether we can write to disk.  If so return 1, if not return 0.

   openw, olun, filename, /get_lun, /delete, error=err
   if (keyword_set(olun)) then close, olun
   if (keyword_set(err)) then begin
      splog, 'Unable to write file to disk: ' + filename
      return, 0B
   endif

   return, 1B
end
;------------------------------------------------------------------------------
function first_readascii, filename, range=range1, columns=columns

   ;----------
   ; Now look for an un-compressed or ASCII version.

   thisfile = (findfile(filename))[0]
   if (NOT keyword_set(thisfile)) then begin
      splog, 'WARNING: FIRST data file not found ' + filename
      splog, 'You may download these data from:'
      splog, '  http://sundog.stsci.edu/first/catalogs/'
      return, 0
   endif

   ;----------
   ; Read a single FIRST data file into a byte array

;  NUMLINES fails if the index file is compressed
;   expectlines = 811117L
   if (keyword_set(range1)) then range = range1 $
    else range = [0L, numlines(thisfile)-1]
;   if (keyword_set(range1)) then range = range1 $
;    else range = [0L, expectlines-1]
   nline = range[1] - range[0] + 1
   barr = bytarr(107,nline)
   tmpString = ''
   get_lun, ilun
   openr, ilun, thisfile
   point_lun, ilun, 108L * range[0]
   for iline=0L, nline-1 do begin
      readf, ilun, tmpString
      barr[0:strlen(tmpString)-1,iline] = byte( tmpString )
   endfor
   close, ilun
   free_lun, ilun

   ; Create data structure and allocate memory for data.
   pTemp = $
   { fobj, $
     first_ra        : 0.d0          , $  ; RA, equinox=2000 [deg]
     first_dec       : 0.d0          , $  ; DEC, equinox=2000 [deg]
     first_warning   : ' '           , $  ; 'W' if sidelobe warning
     first_fint      : 0.0           , $  ;
     first_fpeak     : 0.0           , $  ;
     first_rms       : 0.0           , $  ;
     first_maj       : 0.0           , $  ;
     first_min       : 0.0           , $  ;
     first_pa        : 0.0           , $  ;
     first_fmaj      : 0.0           , $  ;
     first_fmin      : 0.0           , $  ;
     first_fpa       : 0.0           , $  ;
     first_fieldname : '            '  $  ;
   }
   bigdat = replicate( {fobj}, nline )

   ; Copy byte data into data structure
   bigdat.first_ra       = first_convert(barr[  0:  1,*], 'I') * 15.D $
                         + first_convert(barr[  3:  5,*], 'I') * 15.D/60.D $
                         + first_convert(barr[  6: 11,*], 'F') * 15.D/3600.D
   bigdat.first_dec      = (1 - 2*(first_convert(barr[13,*],'A') EQ '-')) * $
                         ( first_convert(barr[ 14: 15,*], 'I') $
                         + first_convert(barr[ 17: 18,*], 'I') /60.D $
                         + first_convert(barr[ 20: 24,*], 'F') /3600.D )
   bigdat.first_warning  = first_convert(barr[ 26: 26,*], 'A')
   bigdat.first_fpeak    = first_convert(barr[ 28: 35,*], 'F')
   bigdat.first_fint     = first_convert(barr[ 37: 45,*], 'F')
   bigdat.first_rms      = first_convert(barr[ 47: 53,*], 'F')
   bigdat.first_maj      = first_convert(barr[ 55: 60,*], 'F')
   bigdat.first_min      = first_convert(barr[ 62: 67,*], 'F')
   bigdat.first_pa       = first_convert(barr[ 69: 73,*], 'F')
   bigdat.first_fmaj     = first_convert(barr[ 75: 80,*], 'F')
   bigdat.first_fmin     = first_convert(barr[ 82: 87,*], 'F')
   bigdat.first_fpa      = first_convert(barr[ 89: 93,*], 'F')
   bigdat.first_fieldname= first_convert(barr[ 95:106,*], 'A')

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
pro first_rdindex, filename, first_dir=first_dir

   common com_first_read, first_index

   ;----------
   ; Find the FIRST index file if it exists.

   indfile = filepath('firstindex.fits', root_dir=first_dir)
   thisfile = (findfile(indfile))[0]

   if (keyword_set(thisfile)) then begin
      first_index = mrdfits(thisfile, 1)
      return
   endif

   ;----------
   ; Generate the index file, with one index per 1000 lines.

   splog, 'Building index file...'
   totline = numlines(filename)
   nper = 1000L
   nindex = long((totline-1) / nper) + 1

   ; Create data structure and allocate memory for data.
   first_index = replicate( $
   { firstind, $
     istart    : 0L , $
     ra        : 0.d  $
   }, nindex)

   ; Loop through the actual data file every NPER line,
   ; and assemble the index list.

   for i=0L, nindex-2 do begin
      iline = i * nper
      onedat = first_readascii(filename, range=[iline,iline])
      first_index[i].istart = iline
      first_index[i].ra = onedat.first_ra
   endfor
   first_index[0].ra = -1.D
   first_index[nindex-1].istart = totline
   first_index[nindex-1].ra = 361.D

   ;----------
   ; If we can write this index file to disk, do so.

   if (first_testwrite(indfile)) then $
    mwrfits, first_index, indfile, /create

   return
end
;------------------------------------------------------------------------------
; If a FITS data file exists, return the name of that file.
; Otherwise, generate the FITS file if it is possible to write to disk,
; and return the name of that newly-written FITS file.
; Otherwise, return null string.

function first_readfits, filename

   ;----------
   ; Return the name of the FITS file (or compressed file) if it exists.

   thisfile = (findfile(filename+'.fits*'))[0]
   if (keyword_set(thisfile)) then return, thisfile

   ;----------
   ; Test whether we can write to disk.  If not, then return a null string.

   fitsfile = filename+'.fits'
   if (first_testwrite(fitsfile) EQ 0) then return, ''

   ;----------
   ; Read the full ASCII file, and write a FITS file.
   ; Return the name of that FITS file.

   splog, 'Generating FITS version of ' + filename
   thisdat = first_readascii(filename)
   if (keyword_set(thisdat)) then begin
      mwrfits, thisdat, fitsfile, /create
   endif else begin
      splog, 'Failure to generate FITS file of ' + filename
      fitsfile = ''
   endelse

   return, fitsfile
end
;------------------------------------------------------------------------------
function first_readfile, filename, range=range, columns=columns

   ;----------
   ; Find a single FIRST data file.  First look for a FITS version.
   ; If it does not exist (and cannot be created), then read the ASCII version.

   thisfile = first_readfits(filename)
   if (keyword_set(thisfile)) then begin
;      splog, 'Reading FITS file ' + thisfile + ' RANGE=', range[0], range[1]
      bigdat = mrdfits(thisfile, 1, range=range, columns=columns, /silent)
   endif else begin
      bigdat = first_readascii(filename, range=range, columns=columns)
   endelse

   return, bigdat
end
;------------------------------------------------------------------------------
function first_read, racen, deccen, radius, $
 node=node, incl=incl, hwidth=hwidth, columns=columns

   common com_first_read, first_index

   first_dir = getenv('FIRST_DIR')
   if (NOT keyword_set(first_dir)) then begin
      splog, 'WARNING: FIRST_DIR environment variable not set.'
      splog, 'This must point to a data directory with files downloade from:'
      splog, '  http://sundog.stsci.edu/first/catalogs/'
      return, 0
   endif

   filename = 'catalog_03apr11.bin'

   ;----------
   ; Read the Tycho index file if it is not already cached in the common block

   if (NOT keyword_set(first_index)) then $
    first_rdindex, filename, first_dir=first_dir
   if (NOT keyword_set(first_index)) then return, 0
   nindex = n_elements(first_index)

   ;----------
   ; Determine which RA range (and row numbers) of the data file to read

   qreadcircle = (n_elements(racen) GT 0 AND n_elements(deccen) GT 0 $
    AND n_elements(radius) GT 0)
   qreadstripe = (n_elements(node) GT 0 AND n_elements(incl) GT 0 $
    AND n_elements(hwidth) GT 0)

   if (qreadcircle) then begin
      decmin = deccen - radius
      decmax = deccen + radius
      if (decmax GE 90) then cosdec = 1.d0 $
       else cosdec = cos(decmax/!radeg)
      ramin = racen - radius / cosdec
      ramax = racen + radius / cosdec
      if (ramax - ramin GE 360 OR decmax GE 90) then begin
         ramin = 0.0
         ramax = 360.0
      endif
      ramin = ramin > 0
      ramax = ramax < 360
   endif else begin
      ramin = 0.0
      ramax = 360.0
   endelse

   i1 = (reverse(where(first_index.ra LE ramin)))[0]
   i2 = (where(first_index.ra GE ramax))[0]
   range = [first_index[i1].istart, first_index[i2].istart-1]

   ;----------
   ; Read the files

   bigdat = first_readfile(filepath(filename, root_dir=first_dir), $
    range=range, columns=columns)

   ;----------
   ; Further trim to the requested circle on the sky

   if (qreadcircle OR qreadstripe) then begin
      ntot = n_elements(bigdat)
      qkeep = bytarr(ntot) + 1B
      if (qreadcircle) then begin
         adiff = djs_diff_angle(bigdat.first_ra, bigdat.first_dec, $
          racen, deccen)
         qkeep = qkeep AND (adiff LE radius)
      endif
      if (qreadstripe) then begin
         radec_to_munu, bigdat.first_ra, bigdat.first_dec, $
          node=node, incl=incl, mu, nu
         qkeep = qkeep AND (abs(nu) LE hwidth)
      endif
      ikeep = where(qkeep, nkeep)
      if (nkeep EQ 0) then bigdat = 0 $
       else bigdat = bigdat[ikeep]
   endif

   return, bigdat
end
;------------------------------------------------------------------------------
