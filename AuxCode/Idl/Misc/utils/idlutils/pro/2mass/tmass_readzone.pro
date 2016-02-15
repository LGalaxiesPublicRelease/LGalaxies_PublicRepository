;+------------------------------------------------------------------------  
; NAME:
;       tmass_readzone
;+------------------------------------------------------------------------  
; PURPOSE:
;       Read 2MASS FITS files given RA range out of one deczone.
;+------------------------------------------------------------------------  
; INPUTS:
;   fitspath  - path to catalogue files (.fits and .acc)
;   zone      - zone number (float, 1/10 degrees)
;   ra0,ra1   - ra limits [deg]
;   prefix    - filename prefix: '2mass' for 2MASS
;
;+------------------------------------------------------------------------  
; OUTPUTS:
;   result    - 2mass data structure defined by the FITS files
;
;+------------------------------------------------------------------------  
; COMMENTS:
;   uses point_lun to skip to requested part of file.  Very fast. 
;
;   Requests are padded by 1/10 the interpolation grid spacing.  This
;     padding is trimmed unless that would yield a null result. 
;
;   Warning - this routine interpolates file index positions
;             and works only if the star distribution is approximately
;             uniform within each interpolation grid patch (which it is).  
;+------------------------------------------------------------------------
; MODIFICATION HISTORY:
;  2003-Jul-14 Taken from usno_read by D. Finkbeiner, Princeton
;-
PRO tmass_readzone, fitspath, zone, ra0, ra1, prefix, result
  
; error trapping
  IF (ra0 LT 0.0) OR (ra0 GT 360.0) OR (ra1 LT 0.0) OR (ra1 GT 360.0) THEN $
    message, 'RA out of range'

  IF (ra0 GE ra1) THEN message, 'ra0 >= ra1'

; pad RA range
  pad = 3.75/10. ; one tenth the index interpolation grid spacing
  raread = [ra0-pad, ra1+pad]
  raread = ((raread) < 360.0) > 0.0

; read .acc file (1-indexed)
  zstr = prefix+string(zone, format='(I4.4)')

  accfile  = djs_filepath(zstr+'.acc', root_dir=fitspath)
  fitsfile = djs_filepath(zstr+'.fits', root_dir=fitspath)
  IF file_test(accfile) NE 1 THEN message, 'cannot find file ' + accfile
  grid = dblarr(3, 96)
  openr, rlun, accfile, /get_lun
  readf, rlun, grid
  free_lun, rlun
  ra_acc = reform(grid[0, *])
  ind    = reform(grid[1, *])
  n      = reform(grid[2, *])
  ntag = n_elements(ra_acc)

; pad acc arrays
  ra = [ra_acc, 24.] *15. ; convert to degrees
  indmax = ind[ntag-1]+n[ntag-1]  ; one-indexed (see below)
  ind = [ind, indmax ]

; get (zero-indexed) offsets in .fits file
  indrange = long(interpol(ind, ra, raread))-1L  ; zero-indexed
  if (indrange[1] LT 0) then stop ; this indicates corrupted .acc file

; read .fits file
  data = mrdfits(fitsfile, 1, range=indrange, /silent)

; trim unwanted RA stars
  
  racat = data.tmass_ra
  good = where((racat LE ra1) AND (racat GE ra0), ct)
  IF ct GT 0 THEN BEGIN 
     result = data[good]
  ENDIF ELSE BEGIN 
     result = data  ; keep padding
  ENDELSE 

  return
end
