;+
; NAME:
;   tmass_reformat_psc
;
; PURPOSE:
;   reformat 2MASS PSC from ASCII to FITS files
;
; CALLING SEQUENCE:
;   tmass_reformat_psc, inpath, outpath
;
; INPUTS:
;   inpath  - path to the gzipped psc files downloaded from IPAC
;   outpath - path to temp dir to write FITS files, one per ASCII file.
;
; COMMENTS:
;   convert the gzipped ASCII files to FITS files. 
;   "tmass_" prepended to tag names for merge with SDSS data. 
;
;   Data Revision 5: April 13, 2003: R. Stiening,
;    Modified for ftp version: April 30, 2003: R. Cutri
;
; REVISION HISTORY:
;   2003-Jul-17  Written by D. Finkbeiner & D. Schlegel, Princeton
;----------------------------------------------------------------------

pro tmass_reformat_psc, inpath, outpath

  spawn, 'mkdir -p '+outpath

  flist = findfile(inpath+'psc*.gz', count=fct)
  infile = strmid(flist, strlen(inpath))

  for i=0L, fct-1 do begin 
     fname = repstr(infile[i], '.gz', '.fits')
     outfile = concat_dir(outpath, fname)
     if file_test(outfile) then begin 
        print, outfile, ' already exists - skipping.'
     endif else begin 
        print, 'reading: ', inpath+infile[i]
        print, 'writing: ', outfile
        spawn, 'touch '+outfile
        tmass_ascii2fits, inpath+infile[i], outfile
     endelse 
  endfor 

  return
end

