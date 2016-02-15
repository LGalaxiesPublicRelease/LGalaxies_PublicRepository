;+
; NAME:
;   beatdrive
;
; PURPOSE:
;   test a new disk partition by writing it full of data and verifying
;
; CALLING SEQUENCE:
;   beatdrive, path, fits=fits, size=size
;
; INPUTS:
;   path   - path to beat on (will make with mkdir -p)
;
; KEYWORDS:
;   fits   - set to write with writefits instead of writeu
;   size   - multiply SDSS image size (5574656 bytes) by this factor
;   
; OUTPUTS:
;   lots of files!
;
; RESTRICTIONS:
;   
; EXAMPLES:
;   beatdrive, '/scr1/test'
;   
; COMMENTS:
;   data file size is 
;   5574656 Bytes = 5444 kB = 5.3164 MB = 0.00519180 GB
;   
; REVISION HISTORY:
;   BOT         Written by Douglas Finkbeiner, Princeton
;   10-Dec-2005 checked into idlutils/pro/misc
;
;----------------------------------------------------------------------
pro beatdrive, path, fits=fits, size=size

; -------- close all files (just in case)
  close, /all

; -------- see if we can make "path"
  file_mkdir, path
  if file_test(path, /dir) EQ 0 then begin 
     print, 'You do not have permission to create ', path
     return
  endif 

; -------- check inputs
  if NOT keyword_set(path) then stop
  if not keyword_set(size) then size = 1
  nperdir = 10000 ; number of files per directory

; -------- pseudo-random image to write (always the same!)
  a = fix(randomu(1, 2048, 1361*size)*10000)
  nmeg = n_elements(a)*2/1048576.

; -------- determine available disk space
  spawn, ['df', '-m', path], res, /noshell
  megfree = strmid(res[1], 40, 10)

; -------- number of files to write (but less than 1 million!)
  nfile = (long(megfree)/5.3164062d -1000) < 1000000L

  if nfile LE 0 then begin 
     print, 'there is not enough space in '+path
     print, 'to do the test'
     return
  endif

  print, 'Writing', nfile, ' files...'

  t1 = systime(1)
  st = t1
  dt = dblarr(nfile)

; -------- loop over files
  for i=0L, nfile-1 do begin 
     t2 = st
     if (i mod nperdir) eq 0 then begin 
        newdir = concat_dir(path, string(i/nperdir, format='(I3.3)'))
        file_mkdir, newdir
        print & print, 'Making new directory: ', newdir
     endif 
     fname = string(i/nperdir, '/idR-00junk-c9-', i mod nperdir, '.fits', $
                    format='(I3.3,A,I5.5,A)')
     if keyword_set(fits) then begin 
        writefits, concat_dir(path, fname), a
     endif else begin 
        openw, wlun, concat_dir(path, fname), /get_lun
        writeu, wlun, a
        free_lun, wlun
     endelse 
     st = systime(1)
     av = (st-t1)/(i+1.d)
     if (i mod 10) eq 0 then print, i, '  avg.', nmeg/(st-t2), av, st-t2
     dt[i] = (st-t2)
  endfor 

; -------- write out timing test
  outname = 'drivetest_'+strmid(getenv('HOST'), 0, 5)+stregex(path,'[0-9]', /extract) $
    +'.fits'
  writefits,  concat_dir(path, outname), dt

; -------- Now verify with beatdrive_check
  beatdrive_check, path

  return
end
