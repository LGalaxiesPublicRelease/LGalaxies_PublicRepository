;+
; NAME:
;   beatdrive_check
;
; PURPOSE:
;   test a new disk partition by reading output of beatdrive.pro
;
; CALLING SEQUENCE:
;   beatdrive_check, path, fits=fits, nocompare=nocompare
;
; INPUTS:
;   path      - path given to beatdrive.pro
;
; KEYWORDS:
;   fits      - set to read with readfits [not implemented]
;   nocompare - turn off comparisons for speed.
;
; EXAMPLES:
;   beatdrive_check, '/scr1/test'
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
pro beatdrive_check, path, fits=fits, nocompare=nocompare

  tbegin = systime(1)
  filesperchunk=1000

  a = indgen(2048, 1361)
  ref = fix(randomu(1, 2048, 1361)*10000)
  nmeg = n_elements(a)*2/1048576.

  if keyword_set(fits) then message, 'sorry, not implemented'

  dirlist = file_search(concat_dir(path, '???'), /test_dir, count=ndir)

  for idir=0, ndir-1 do begin 
     dir = dirlist[idir]
     print, 'Working on ', dir
     flist = file_search(concat_dir(dir, '/*.fits'), count=fct)
     t1 = systime(1) 
     st = t1

     nchunk=long(fct/filesperchunk)
     for ichunk=0, nchunk-1 do begin 
        t2 = st
        for i=ichunk*filesperchunk,(ichunk+1)*filesperchunk-1 do begin 
            fname = flist[i]
            openr, rlun, fname, /get_lun
            readu, rlun, a
            free_lun, rlun
            if NOT keyword_set(nocompare) then begin 
                if array_equal(a, ref) eq 0 then begin 
                    print, 'In file ', fname
                    message, 'Array not equal to expected reference array.'
                endif
            endif
        endfor
        st = systime(1)
        av = (st-t1)/(i+1.d)
        print, i, '  avg.', nmeg*filesperchunk/(st-t2), av, (st-t2)/filesperchunk
     endfor 
     print, 'files read:', fct
     print, 'In directory: ', dir
  endfor 
  print, 'Ndir: ', strcompress(ndir)
  print, systime()
  Telapse = systime(1)-tbegin
  print, 'Elapsed time', Telapse, ' seconds'
  print, 'Elapsed time', Telapse/3600, ' hours'
  print, 'Finished.'

  return
end
