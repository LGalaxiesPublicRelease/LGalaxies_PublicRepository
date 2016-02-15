;+
; NAME:
;   sdss_flagval
;
; PURPOSE:
;   Return bitmask values corresponding to labels.
;
; CALLING SEQUENCE:
;   value = sdss_flagval(flagprefix, label)
;
; INPUTS:
;   flagprefix - Flag name (scalar string).  The following are supported:
;                SPPIXMASK, TARGET, TTARGET.
;   label      - String name(s) corresponding to each non-zero bit in FLAGVALUE.
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   value      - Signed long with any number of its bits set.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This function is the inverse of SDSS_FLAGNAME().
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   splog
;   yanny_free
;   yanny_read
;
; DATA FILES:
;   $IDLUTILS_DIR/data/sdss/sdssMaskbits.par
;
; REVISION HISTORY:
;   02-Apr-2002 Written by D. Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function sdss_flagval, flagprefix, inlabel

   ; Declare a common block so that the mask names are remembered between calls.
   common com_maskbits, maskbits

   if (n_params() NE 2 OR n_elements(flagprefix) NE 1) then begin
      print, 'Syntax - value = sdss_flagval(flagprefix, label)'
      return, ''
   endif

   ;----------
   ; Read the parameter file the 1st time this function is called.
   ; (After that, store this info in a common block.)

   if (NOT keyword_set(maskbits)) then begin
      maskfile = filepath('sdssMaskbits.par', $
       root_dir=getenv('IDLUTILS_DIR'), subdirectory='data/sdss')
      if (NOT keyword_set(maskfile)) then $
       message, 'File with mask bits not found'
      yanny_read, maskfile, pdat
      maskbits = *pdat[0]
      yanny_free, pdat
   endif

   ;----------
   ; Generate a list of all non-blank labels as a string array

   flagvalue = 0L

   alllabel = strsplit(inlabel[0], /extract)
   for i=1, n_elements(inlabel)-1 do $
    alllabel = [alllabel, strsplit(inlabel[i], /extract)]
   ilabel = where(alllabel NE '', nlabel)
   if (nlabel EQ 0) then return, flagvalue
   alllabel = alllabel[ilabel]

   ;----------
   ; Find the match for each label, and add its value to the output

   for ilabel=0, nlabel-1 do begin
      imatch = where(strupcase(flagprefix[0]) EQ maskbits.flag $
       AND strupcase(alllabel[ilabel]) EQ strupcase(maskbits.label), ct)
      if (ct NE 1) then $
       message, 'ABORT: Unknown bit label ' + strupcase(alllabel[ilabel]) $
        + ' for flag ' + strupcase(flagprefix)

      flagvalue = flagvalue + 2L^(maskbits[imatch[0]].bit)
   endfor

   return, flagvalue
end
;------------------------------------------------------------------------------
