pro ftab_help,file_or_fcb,EXTEN_NO = exten_no, TEXTOUT = textout
;+
; NAME:
;       FTAB_HELP
; PURPOSE:
;       Describe the columns of a FITS binary or ASCII table extension.
;
; CALLING SEQUENCE:
;       FTAB_HELP, filename, [ EXTEN_No = , TEXTOUT= ]
;               or
;       FTAB_HELP, fcb, [EXTEN_No=, TEXTOUT= ]
;
; INPUTS:
;       filename - scalar string giving name of the FITS file.  
;       fcb - FITS control block returned by a previous call to FITS_OPEN
;
; OPTIONAL KEYWORD INPUTS:
;       EXTEN_NO - integer scalar specifying which FITS extension to read.
;               Default is to display the first FITS extension.
;       TEXTOUT - scalar number (0-7) or string (file name) determining
;               output device (see TEXTOPEN).  Default is TEXTOUT=1, output 
;               to the user's terminal    
;
; EXAMPLE:
;       Describe the columns in the second extension of a FITS file spec.fits
;       and write the results to a file 'spec2.lis'
;
;       IDL> ftab_help,'spec.fits',exten=2,t='spec2.lis'
;
; SYSTEM VARIABLES:
;       Uses the non-standard system variables !TEXTOUT and !TEXTUNIT
;       which must be defined (e.g. with ASTROLIB) before compilation
; PROCEDURES USED:
;       FITS_READ, FITS_CLOSE, FITS_OPEN, FTHELP, TBHELP, TEXTOPEN, TEXTCLOSE
; HISTORY:
;       version 1  W. Landsman    August 1997
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Corrected documentation W. Landsman   September 1997
;       Don't call fits_close if fcb supplied W. Landsman May 2001 
;-
;----------------------------------------------------------------------
 if N_params() LT 1 then begin
        print,'Syntax - FTAB_HELP, fcb_or_filename, [EXTEN_NO=, TEXTOUT= ]'
        return
 endif

 if not keyword_set(exten_no) then exten_no = 1

 sz = size(file_or_fcb)                                                    
 if sz[sz[0]+1] NE 8 then fits_open,file_or_fcb,fcb else fcb=file_or_fcb
 if fcb.nextend EQ 0 then begin 
          message,'File contains no Table extensions',/INF
          if sz[sz[0]+1] NE 8 then fits_close,fcb else $
                      file_or_fcb.last_extension = exten_no
          return
  endif

 fits_read,fcb, dummy, htab, /header_only,/no_pdu, exten_no=exten_no
 if sz[sz[0]+1] NE 8 then fits_close,fcb else $
         file_or_fcb.last_extension = exten_no
 ext_type = fcb.xtension[exten_no]

 case ext_type of
 'A3DTABLE': binary = 1b
 'BINTABLE': binary = 1b
 'TABLE': binary = 0b
 else: message,'ERROR - Extension type of ' + $
                ext_type + ' is not a FITS table format'
 endcase

 textopen,'ftab_help',textout=textout
 printf,!TEXTUNIT, 'FITS file: ' + fcb.filename 
 printf,!TEXTUNIT, 'Extension No: ' + strtrim(exten_no,2)

 if binary then tbhelp, htab, TEXTOUT = 5 $
           else fthelp, htab, TEXTOUT = 5
 textclose, textout=textout
 return
 end
