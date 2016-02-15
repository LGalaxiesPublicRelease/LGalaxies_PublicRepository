;+
; NAME:
;   gsc_nsky
;
; PURPOSE:
;   Generate IDL structure for GSC catalogue at dec > -30 and m < 18.0
;
; CALLING SEQUENCE:
;   gsc_nsky, [ catall ]
; 
; INPUTS:
;   <data files in $PHOTO_DATA/gsc>
;
; OUTPUTS:
;   <data file (binary FITS table) in $PHOTO_DATA/kocat>
;
; OPTIONAL OUTPUTS:
;   catall  - Structure containing all objects dec > -30 deg
;
; INPUT FILES:
;   $GSC_DIR/n????/*.gsc
;   $GSC_DIR/s????/*.gsc
;
; OUTPUT FILES:
;   $GSC_DIR/gsc-Nsky.fits
;
; EXAMPLES:
;
; COMMENTS:
;   Note: these catalogues are written as F9.5 ASCII FITS tables.  So 
;    there is no reason to read them as double precision(?)
;   For more information see gsc/readme.txt
;   This can take a long time to run. 
;
; PROCEDURES CALLED:
;   djs_angle_match()
;   gsc_read_table()
;   mwrfits
;
; REVISION HISTORY:
;   Written 2002-Mar-27 by D. Finkbeiner
;-
;------------------------------------------------------------------------------
function gsc_dec2hms_arr, decarr

   n = n_elements(decarr)
   s = strarr(n)

   for i=0L, n-1 do begin 
      s[i] = dec2hms(decarr[i])
   endfor

   return, s
end
;------------------------------------------------------------------------------
pro gsc_nsky, catall

   delvarx, catall
   gsc_dir = getenv('GSC_DIR')
   if (NOT keyword_set(gsc_dir)) then $
    message, 'GSC_DIR not set!'
   maglim = [0, 18.0]
   poserrlim = 0.55  ; arcsec

   ;----------
   ; Construct the subdirectory list within $GSC_DIR

   hms = strmid(strcompress(gsc_dec2hms_arr(findgen(12)*7.5), /rem), 0, 4)
   nlist = 'n'+hms
   slist = 's'+hms
   dirlist = [nlist, slist[0:3]] ; ignore some of the south

   ;----------
   ; Loop through each directory, reading all files

   for j=0L, n_elements(dirlist)-1 do begin 
      flist = findfile(filepath('*.gsc', root_dir=gsc_dir, subdir=dirlist[j]))
      print, 'Directory', dirlist[j], '  ', systime()

      cat = gsc_read_table(flist, maglim=maglim, poserrlim=poserrlim)

      catall = keyword_set(catall) ? [catall, cat] : cat
   endfor 

   ;----------
   ; Keep only those stars without neighbors within 4 arcsec.
   ; (These are mostly catalogue duplicates, not real duplicates!

   ntot = djs_angle_match(catall.ra, catall.dec, dtheta=4./3600, $
    mcount=mcount,mdist=mdist)
   indx = where(mcount EQ 0)
   catall = catall[indx]

   ;----------
   ; Write the output file

   mwrfits, catall, filepath('gsc-Nsky.fits', root_dir=gsc_dir), /create

   return
end
;------------------------------------------------------------------------------
