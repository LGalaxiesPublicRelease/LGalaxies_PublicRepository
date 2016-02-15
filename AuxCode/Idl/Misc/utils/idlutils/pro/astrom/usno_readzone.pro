;+
;+------------------------------------------------------------------------  
; NAME:
;       usno_readzone
;+------------------------------------------------------------------------  
; PURPOSE:
;       Read given RA range out of one deczone. (USNO-A or B)
;+------------------------------------------------------------------------  
; INPUTS:
;   catpath   - path to catalogue files (.cat and .acc)
;   zone      - zone number (float, 1/10 degrees)
;   ra0,ra1   - ra limits [deg]
;   rec_len   - record length for each object [bytes] (read as ulongs)
;   prefix    - filename prefix: 'zone' for USNOA, 'b' for USNOB
;
; KEYWORDS:
;   swap_if... - USNOA is written little endian, B is big endian(!)
;+------------------------------------------------------------------------  
; OUTPUTS:
;   data      - float(rec_len/4,N) array of results.  
;                 data[0,*] = RA (in .01 arcsec)
;                 data[1,*] = (dec+90) (in .01 arcsec)
;       USNO-A:
;                 data[2,*] = magnitudes packed in 32-bit int
;       USNO-B:   
;                 data[2:19,*] = all kinds of stuff
;+------------------------------------------------------------------------  
; COMMENTS:
;   uses point_lun to skip to requested part of file.  Very fast. 
;
;   Requests are padded by 1/10 the interpolation grid spacing.  This
;     padding is trimmed unless that would yield a null result. 
;
;   Warning - this routine interpolates file index positions
;             and works only if the star distribution is approximately
;             uniform (which it is).  
;+------------------------------------------------------------------------
; MODIFICATION HISTORY:
;  2002-Nov-25 Taken from usno_read by D. Finkbeiner, Princeton
;               Also works for USNO-B1.0 now. 
;-
PRO usno_readzone, catpath, zone, ra0, ra1, rec_len, prefix, result, $
        swap_if_little_endian=swap_if_little_endian, $
        swap_if_big_endian=swap_if_big_endian
  
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

  accfile = djs_filepath(zstr+'.acc', root_dir=catpath)
  catfile = djs_filepath(zstr+'.cat', root_dir=catpath)
  flist = findfile(accfile, count=ct)
  IF ct NE 1 THEN message, 'cannot find file ' + accfile
  grid = dblarr(3, 96)
  openr, rlun, accfile, /get_lun
  readf, rlun, grid
  free_lun, rlun
  ra_acc = reform(grid[0, *])
  ind    = reform(grid[1, *])
  n      = reform(grid[2, *])
  ntag = n_elements(ra_acc)

; pad arrays
  ra = [ra_acc, 24.] *15. ; convert to degrees
  indmax = ind[ntag-1]+n[ntag-1]-1
  ind = [ind, indmax ]

; get (zero-indexed) offsets in .cat file
  indrange = long(interpol(ind, ra, raread))-1L
  ind0 = indrange[0]
  ind1 = indrange[1]

; read .cat file
  openr, readlun, catfile, /get_lun, $
    swap_if_little_endian=swap_if_little_endian, $
    swap_if_big_endian=swap_if_big_endian
  nstars = (ind1-ind0)+1
  data = ulonarr(rec_len/4, nstars)
  point_lun, readlun, ind0*rec_len
  readu, readlun, data
  free_lun, readlun

; trim unwanted RA stars
  racat  = transpose(data[0, *]) /3.6d5
  good = where((racat LE ra1) AND (racat GE ra0), ct)
  IF ct GT 0 THEN BEGIN 
     result = data[*, good]
  ENDIF ELSE BEGIN 
     result = data  ; keep padding
  ENDELSE 

  return
end
