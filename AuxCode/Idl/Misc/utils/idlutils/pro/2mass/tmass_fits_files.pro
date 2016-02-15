;+
; NAME:
;   tmass_fits_files
;
; PURPOSE:
;   Convert 2MASS gzipped ASCII files to decsliced FITS tables
;
; CALLING SEQUENCE:
;   tmass_fits_files
;
; INPUTS:
;   <hardwired paths>
;
; OPTIONAL INPUTS:
;   
; KEYWORDS:
;   
; OUTPUTS:
;   
; OPTIONAL OUTPUTS:
;   
; RESTRICTIONS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   The first step, tmass_reformat_psc, can be run on several boxes
;    in parallel. 
;
; REVISION HISTORY:
;   2003-Jul-24  Written by Douglas Finkbeiner, Princeton
;----------------------------------------------------------------------

pro tmass_fits_files

; -------- convert ASCII to FITS
  ascii_dir = '/peyton/scr/photo66/sdssdata/2mass/allsky/'
  temp_dir  = '/peyton/scr/photo63/sdssdata/2mass/allsky/'
  tmass_reformat_psc, ascii_dir, temp_dir

; -------- write declination slices
  twomass_dir = concat_dir(getenv('TWOMASS_DIR'), 'slice')
  if twomass_dir eq 'slice' then begin 
     print, 'You must set TWOMASS_DIR'
     stop
  endif 

  tmass_decslice, temp_dir, twomass_dir



  return
end

