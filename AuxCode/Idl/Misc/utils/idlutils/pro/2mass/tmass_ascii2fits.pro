;+
; NAME:
;   tmass_ascii2fits
;
; PURPOSE:
;   Convert 2MASS ascii files to FITS tables
;
; CALLING SEQUENCE:
;   tmass_ascii2fits, infile, outfile
;
; INPUTS:
;   infile  - input file name
;   outfile - output file name
;
; KEYWORDS:
;   
; OUTPUTS:
;   file named <outfile>
;
; RESTRICTIONS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   Recognizes gzip compression automatically. 
;
;   I'm not keeping all the fields - just the useful ones. 
;   Might want to keep everything in the future. 
;
; REVISION HISTORY:
;   2003-Jun-26  Written by Douglas Finkbeiner, Princeton
;
;----------------------------------------------------------------------
pro tmass_ascii2fits, infile, outfile

; -------- check for input file
  if file_test(infile) eq 0 then begin 
     print, 'file ', infile, ' not found'
     return
  endif

  gz = strmid(infile, strlen(infile)-3, 3) eq '.gz'
  t0 = systime(1)
  maxline = 50000
  print
  print, 'Scanning ', infile
  if gz then begin 
     spawn, 'gunzip < '+infile+' | wc -l ', result
     openr, rlun, infile, /get_lun, /compress
  endif else begin 
     spawn, 'wc -l '+infile, result
     openr, rlun, infile, /get_lun
  endelse 

  nrec = long((strsplit(result, ' ', /extract))[0])
  print, nrec, ' records found'
  fptr = 0L

  a = tmass_struc(nrec)

  while fptr lt nrec-1 do begin
     nline = (nrec-fptr) < maxline
     str = strarr(nline)
     readf, rlun, str
     print, fptr, nline, systime(1)-t0
     for j=0L, n_elements(str)-1 do begin 
        b = a[0]
        field = strsplit(str[j], '|', /extract)
        field = repstr(field, '\N', '0')
        magerr = float( field[[7,11,15]] )
        b.tmass_ra      = field[0]
        b.tmass_dec     = field[1]
        b.tmass_err_maj = field[2]
        b.tmass_err_min = field[3]
        b.tmass_err_ang = field[4]
        b.tmass_j       = field[6]
        b.tmass_j_ivar  = (magerr[0] LE 0) ? 0 : 1./magerr[0]^2
        b.tmass_h       = field[10]
        b.tmass_h_ivar  = (magerr[1] LE 0) ? 0 : 1./magerr[1]^2
        b.tmass_k       = field[14]
        b.tmass_k_ivar  = (magerr[2] LE 0) ? 0 : 1./magerr[2]^2
        b.tmass_ph_qual = field[18]
        b.tmass_rd_flg  = field[19]
        b.tmass_bl_flg  = field[20]
        b.tmass_cc_flg  = field[21]
        detstr = field[22]
        b.tmass_ndetect = $
         fix([strmid(detstr,0,1), strmid(detstr,2,1), strmid(detstr,4,1)])
        b.tmass_nobserve = $
         fix([strmid(detstr,1,1), strmid(detstr,3,1), strmid(detstr,5,1)])
        b.tmass_gal_contam = field[26]
        b.tmass_mp_flg = field[27]
        b.tmass_pts_key = field[28]
        b.tmass_hemis = field[29]
        b.tmass_jdate = field[35]
        b.tmass_dup_src = field[48]
        b.tmass_use_src = field[49]
        a[fptr+j] = b

     endfor 

     fptr = fptr+nline
  endwhile
  
  free_lun, rlun

  print, 'Writing ', outfile
  mwrfits_chunks, a, outfile, chunksize=50000, /create

  return
end

