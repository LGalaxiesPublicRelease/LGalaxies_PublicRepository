;+
; NAME:
;   tmass_decslice
;
; PURPOSE:
;   rewrite 2MASS data in 0.1 deg declination slices for fast seeks
;
; CALLING SEQUENCE:
;   
; INPUTS:
;   inpath  - path of input files
;   outpath - path of output files
;
; COMMENTS:
;   called by tmass_ingest
;
; REVISION HISTORY:
;   2003-July-12  Written by Douglas Finkbeiner, Princeton
;   2005-June-25  Fix bug that corrupts .acc file if zero stars in any
;           RA step.  - DPF
;
;----------------------------------------------------------------------

pro tmass_writeslice, ind, a, outpath

; -------- sort list on RA
  sind = sort(a.tmass_ra)
  a = a[sind]
  delvarx, sind

; -------- write FITS file
  fname = string('2mass-', ind, format='(A,I4.4)')
  path = concat_dir(outpath, string(ind/10, format='(I3.3)'))+'/'
  spawn, 'mkdir -p '+ path
  print, 'Writing '+path+fname+'.fits'
  mwrfits, a, path+fname+'.fits', /create

; -------- write index (acc) file
  rah = a.tmass_ra/15
  step = 0.25d
  openw, wlun, path+fname+'.acc', /get_lun
  thisind = 1
  for j=0, round(24/step)-1 do begin 
     jnd = where((rah ge (j*step)) and (rah lt ((j+1)*step)), ngood)
     printf, wlun, j*step, thisind, ngood, format='(F5.2,I12,I12)'
     thisind += ngood
  endfor 
  free_lun, wlun

  return
end



pro tmass_decslice, inpath, outpath

  flist = findfile(inpath+'psc_*.fits')
  fileind = 0
  print, 'Reading ', flist[fileind]
  a = mrdfits(flist[fileind], 1)
  ind = ceil((a.tmass_dec+90)*10)-1  ; this is how they did it in the south
  if min(ind) lt 0 then stop
  i = min(ind, max=mx)
  while (i lt 1800) do begin 
     w = where(ind eq i, ngood)
     if ngood le 1 then message, 'no objects!'
     aw = a[w]
     tmass_writeslice, i, aw, outpath
     i = i+1
; -------- check whether we need to read in next file
     if i eq mx then begin 
        w = where(ind eq mx)
        aw = a[w]
        delvarx, a
        fileind = fileind+1
        if fileind lt n_elements(flist) then begin 
           print, 'Reading ', flist[fileind]
           a = [aw, mrdfits(flist[fileind], 1)]
           ind = ceil((a.tmass_dec+90)*10)-1
           if min(a.tmass_dec) ge 0 then ind = floor((a.tmass_dec+90)*10)
           mn = min(ind, max=mx)
           if i ne mn then stop ; assert i=mn
        endif else a = aw
     endif 

  endwhile

  return
end
