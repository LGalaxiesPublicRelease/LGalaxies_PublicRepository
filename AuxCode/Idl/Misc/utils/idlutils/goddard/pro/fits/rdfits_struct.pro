pro rdfits_struct, filename, struct,SILENT = silent, HEADER_ONLY = header_only 
;+
; NAME:
;      RDFITS_STRUCT
; PURPOSE:
;      Read an entire FITS file (all extensions) into a single IDL structure. 
; EXPLANATION:
;      Each header, image or table array is placed in a separate structure 
;      tag.
;
; CALLING SEQUENCE:
;      RDFITS_STRUCT, filename, struct, /SILENT, /HEADER_ONLY ]
;
; INPUT:
;      FILENAME = Scalar string giving the name of the FITS file.  
;                 One can also specify a gzip (.gz) compressed file 
;
; OPTIONAL KEYWORD: 
;      /HEADER_ONLY - If set, then only the FITS headers (and not the data)
;                are read into the structure.
;      /SILENT - Set this keyword to suppress informational displays at the
;               terminal.
; OUTPUT:
;      struct = structure into which FITS data is read.   The primary header
;             and image are placed into tag names HDR0 and IM0.   The ith
;             extension is placed into the tag names HDRi, and either TABi
;             (if it is a binary or ASCII table) or IMi (if it is an image
;             extension)
;
;             If /HEADER_ONLY is set, then struct will contain tags HDR0, HDR1
;             ....HDRn containing all the headers of a FITS file with n 
;             extensions
; PROCEDURES USED:
;       FITS_OPEN, FITS_READ, FITS_CLOSE
;
; METHOD:
;       The file is opened with FITS_OPEN which return information on the 
;       number and type of each extension.    The CREATE_STRUCT() function
;       is used iteratively, with FITS_READ calls to build the final structure.
;
; EXAMPLE:
;       Read the FITS file 'm33.fits' into an IDL structure, st
;
;       IDL> rdfits_struct, 'm33.fits', st
;       IDL> help, /str, st                   ;Display info about the structure
;
; RESTRICTIONS:
;       Does not handle random groups or variable length binary tables
; MODIFICATION HISTORY:
;       Written K. Venkatakrishna, STX April 1992
;       Code cleaned up a bit  W. Landsman  STX  October 92
;       Modified for MacOS     I.  Freedman  HSTX April 1994
;       Work under Windows 95  W. Landsman   HSTX  January 1996
;       Use anonymous structures, skip extensions without data WBL April 1998
;       Converted to IDL V5.0, W. Landsman, April 1998
;       OS-independent deletion of temporary file  W. Landsman  Jan 1999
;       Major rewrite to use FITS_OPEN and CREATE_STRUCT() W. Landsman Sep 2002
;       Added /HEADER_ONLY keyword   W. Landsman  October 2003
;       Do not copy primary header into extension headers W. Landsman Dec 2004
;       Do not modify NAXIS when using /HEADER_ONLY W. Landsman Jan 2005
;-

 compile_opt idl2
 if N_Params() LT 2 then begin 
        print,'Syntax - RDFITS_STRUCT, file, struct, [ /SILENT, /HEADER_ONLY ]'
        return
 endif

 fits_open, filename, fcb                ; Get the description of the file
 if not keyword_set(silent) then $
      message,/inf,'Now reading file ' + filename + ' with ' + $
      strtrim(fcb.nextend,2) + ' extensions'

 h_only = keyword_set(header_only)  
 if h_only then begin
     fits_read,fcb,0,h,/header_only,exten_no=0
     struct = {hdr0:h}
 endif else begin
     fits_read,fcb,d,h,exten_no=0
     struct = {hdr0:h,im0:temporary(d)}
 endelse

 if fcb.nextend EQ 0 then begin
      fits_close,fcb 
      return
 endif
 for i=1,fcb.nextend do begin
     if h_only then begin
     fits_read,fcb,0,h,/header_only,/no_pdu
     struct = create_struct(temporary(struct), 'hdr' + strtrim(i,2), $
              temporary(h))
     endif else begin
     fits_read,fcb,d,h,/no_pdu
     if fcb.xtension[i] EQ 'IMAGE' then tag = 'im' + strtrim(i,2) $
                                else tag = 'tab' + strtrim(i,2)
     struct = create_struct(temporary(struct), 'hdr' + strtrim(i,2), $
              temporary(h),tag, temporary(d))
    endelse
 endfor
     
 fits_close,fcb                             
 return
 end
