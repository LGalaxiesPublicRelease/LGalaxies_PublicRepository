pro MODFITS, filename, data, header, EXTEN_NO = exten_no, ERRMSG = errmsg
           
;+
; NAME:
;      MODFITS
; PURPOSE:
;      Modify a FITS file by updating the header and/or data array.  
; EXPLANATION:
;      Since August 2003, the size of the supplied FITS header or data array
;      does not need to match the size of the existing header or data array.
;
; CALLING SEQUENCE:
;      MODFITS, Filename_or_fcb, Data, [ Header, EXTEN_NO =, ERRMSG=]
;
; INPUTS:
;      FILENAME/FCB = Scalar string containing either the name of the FITS file  
;                  to be modified, or the IO file control block returned after 
;                  opening the file with FITS_OPEN,/UPDATE.   The explicit
;                  use of FITS_OPEN can save time if many extensions in a 
;                  single file will be updated.
;
;      DATA - data array to be inserted into the FITS file.   Set DATA = 0
;               to leave the data portion of the FITS file unmodified
;
;      HEADER - FITS header (string array) to be updated in the FITS file.
;
; OPTIONAL INPUT KEYWORDS:
;      EXTEN_NO - scalar integer specifying the FITS extension to modified.  For
;               example, specify EXTEN = 1 or /EXTEN to modify the first 
;               FITS extension. 
; OPTIONAL OUTPUT KEYWORD:
;       ERRMSG - If this keyword is supplied, then any error mesasges will be
;               returned to the user in this parameter rather than depending on
;               on the MESSAGE routine in IDL.   If no errors are encountered
;               then a null string is returned.               
;
; EXAMPLES:
;     (1) Modify the value of the DATE keyword in the primary header of a 
;             file TEST.FITS.
;
;              IDL> h = headfits('test.fits')      ;Read primary header
;              IDL> sxaddpar,h,'DATE','2001-03-23' ;Modify value of DATE 
;              IDL> modfits,'test.fits',0,h        ;Update header only
;
;       (2) Replace the values of the primary image array in 'test.fits' with 
;               their absolute values
;
;               IDL> im = readfits('test.fits')    ;Read image array
;               IDL> im = abs(im)                  ;Take absolute values
;               IDL> modfits,'test.fits',im        ;Update image array
;
;       (3) Add some HISTORY records to the FITS header in the first extension
;               of a file 'test.fits'
;       
;               IDL> h = headfits('test.fits',/ext)  ;Read first extension hdr
;               IDL> sxaddhist,['Comment 1','Comment 2'],h
;               IDL> modfits,'test.fits',0,h,/ext    ;Update extension hdr
;
;       (4) Change 'OBSDATE' keyword to 'OBS-DATE' in every extension in a 
;           FITS file.    Explicitly open with FITS_OPEN to save compute time.
;
;               fits_open,'test.fits',io,/update    ;Faster to explicity open
;               for i = 1,nextend do begin          ;Loop over extensions
;                   fits_read,io,0,h,/header_only,exten_no=i,/No_PDU ;Get header     
;                   date= sxpar(h,'OBSDATE')         ;Save keyword value
;                   sxaddpar,h,'OBS-DATE',date,after='OBSDATE' 
;                   sxdelpar,h,'OBSDATE'             ;Delete bad keyword
;                   modfits,io,0,h,exten_no=i        ;Update header
;               endfor
;
;           Note the use of the /No_PDU keyword in the FITS_READ call -- one 
;           does *not* want to append the primary header, if the STScI 
;           inheritance convention is adopted.
;
; NOTES:
;       Uses the BLKSHIFT procedure to shift the contents of the FITS file if 
;       the new data or header differs in size by more than 2880 bytes from the
;       old data or header.
;
;       Also see the procedures FXHMODIFY to add a single FITS keyword to a 
;       header in a FITS files, and FXBGROW to enlarge the size of a binary 
;       table.
; RESTRICTIONS:
;       (1) Cannot be used to modifiy the data in FITS files with random 
;           groups or variable length binary tables.   (The headers in such
;           files *can* be modified.)
;
;       (2) If a data array but no FITS header is supplied, then MODFITS does 
;           not check to make sure that the existing header is consistent with
;           the new data.
;
;       (3) Does not work with compressed files
; PROCEDURES USED:
;       Functions:   IS_IEEE_BIG(), N_BYTES(), SXPAR()
;       Procedures:  BLKSHIFT, CHECK_FITS, FITS_OPEN, FITS_READ, HOST_TO_IEEE
;
; MODIFICATION HISTORY:
;       Written,    Wayne Landsman          December, 1994
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Fixed possible problem when using WRITEU after READU   October 1997
;       New and old sizes need only be the same within multiple of 2880 bytes
;       Added call to IS_IEEE_BIG()     W. Landsman   May 1999
;       Added ERRMSG output keyword     W. Landsman   May 2000
;       Update tests for incompatible sizes   W. Landsman   December 2000
;       Major rewrite to use FITS_OPEN procedures W. Landsman November 2001
;       Add /No_PDU call to FITS_READ call  W. Landsman  June 2002
;       Update CHECKSUM keywords if already present in header, add padding 
;       if new data  size is smaller than old  W.Landsman December 2002
;       Only check XTENSION value if EXTEN_NO > 1   W. Landsman Feb. 2003
;       Correct for unsigned data on little endian machines W. Landsman Apr 2003
;       Major rewrite to allow changing size of data or header W.L. Aug 2003
;       Fixed case where updated header exactly fills boundary W.L. Feb 2004
;       More robust error reporting W.L. Dec 2004
;-
  On_error,2                    ;Return to user

; Check for filename input

   if N_params() LT 1 then begin                
      print,'Syntax - MODFITS, Filename, Data, [ Header, EXTEN_NO=, ERRMSG= ]'
      return
   endif

   if not keyword_set( EXTEN_NO ) then exten_no = 0
   if N_params() LT 2 then Header = 0
   nheader = N_elements(Header)
   ndata = N_elements(data)
   dtype = size(data,/TNAME)
   printerr =  not arg_present(ERRMSG) 
   fcbsupplied = size(filename,/TNAME) EQ 'STRUCT'

   if (nheader GT 1) and (ndata GT 1) then begin
        check_fits, data, header, /FITS, ERRMSG = MESSAGE
        if message NE '' then goto, BAD_EXIT
   endif

; Open file and read header information
         
   if exten_no EQ 0 then begin 
         if nheader GT 0 then $
             if strmid( header[0], 0, 8)  NE 'SIMPLE  ' then begin 
                 message = $
                'Input header does not contain required SIMPLE keyword'
                 goto, BAD_EXIT
             endif
   endif else begin
         if nheader GT 0 then $
             if strmid( header[0], 0, 8)  NE 'XTENSION' then begin 
              message = $
             'ERROR - Input header does not contain required XTENSION keyword'
              goto, BAD_EXIT
              endif
   endelse

; Was a file name or file control block supplied?

   if not fcbsupplied then begin 
       fits_open, filename, io,/update,/No_Abort,message=message
       if message NE '' then GOTO, BAD_EXIT
    endif else begin 
       if filename.open_for_write EQ 0 then begin
             message = 'FITS file is set for READONLY, cannot be updated'
             goto, BAD_EXIT
       endif
       io = filename
   endelse

; Getting starting position of data and header

   unit = io.unit
   start_d = io.start_data[exten_no]
   if exten_no NE io.nextend then begin
        start_h = io.start_header[exten_no+1] 
        nbytes = start_h - start_d
   endif else nbytes = N_BYTES(data)

   FITS_READ,Io,0,oldheader,/header_only, exten=exten_no, /No_PDU, $
       message = message,/no_abort
   if message NE '' then goto, BAD_EXIT
    dochecksum = sxpar(oldheader,'CHECKSUM', Count = N_checksum)
   checksum = N_checksum GT 0  

; Update header, including any CHECKSUM keywords if present 

   if nheader GT 1 then begin
        noldheader = start_d - io.start_header[exten_no]
 
        if dtype EQ 'UINT' then $
              sxaddpar,header,'BZERO',32768,'Data is unsigned integer'
        if dtype EQ 'ULONG' then $
              sxaddpar,header,'BZERO',2147483648,'Data is unsigned long'
        if checksum then begin 
               if Ndata GT 1 then FITS_ADD_CHECKSUM, header, data else $
                FITS_ADD_CHECKSUM, header 
        endif
        newbytes = N_elements(header)*80 
        block = (newbytes-1)/2880 - (Noldheader-1)/2880
        if block NE 0 then begin  
            BLKSHIFT, io.unit, start_d, block*2880L
            start_d = start_d + block*2880L
            if exten_no NE io.nextend then start_h = start_h + block*2880L
        endif
        point_lun, unit, io.start_header[exten_no]      ;Position header start  
        writeu, unit, byte(header)
        remain = newbytes mod 2880
	if remain GT 0 then writeu, unit, replicate( 32b, 2880 - remain)
 
   endif 

   if ndata GT 1 then begin
 
        newbytes = N_BYTES(data)    ;total number of bytes in supplied data
        nblock = (newbytes-1)/2880 - (nbytes-1)/2880
        if nblock NE 0 and exten_no NE io.nextend then $
              BLKSHIFT, io.unit, start_h, nblock*2880L
      
        if nheader EQ 0 then begin
                check_fits,data,oldheader,/FITS,ERRMSG = message
                if message NE '' then goto, BAD_EXIT
        endif
        vms = !VERSION.OS EQ "vms"
        Little_endian = not IS_IEEE_BIG()

        junk = fstat(unit)   ;Need this before changing from READU to WRITEU
        point_lun, unit, start_d
        if dtype EQ 'UINT' then newdata = fix(data - 32768)
        if dtype EQ 'ULONG' then newdata = long(data - 2147483648)
        if (VMS or Little_endian) then begin
             if (dtype NE 'UINT') and (dtype NE 'ULONG') then newdata = data
             HOST_TO_IEEE, newdata
        endif
        if N_elements(newdata) GT 0 then writeu, unit, newdata  else $
                                         writeu, unit ,data
        remain = newbytes mod 2880
	if remain GT 0 then begin
             padnum = 0b
             if exten_no GT 0 then begin 
                 exten = sxpar( oldheader, 'XTENSION')
	         if exten EQ 'TABLE   ' then padnum = 32b
             endif
	     writeu, unit, replicate( padnum, 2880 - remain)
	endif
    endif       

   if not fcbsupplied then FITS_CLOSE,io
   return 

BAD_EXIT:
    if N_elements(io) GT 0 then if not fcbsupplied then fits_close,io
    if printerr then message,'ERROR - ' + message,/CON else errmsg = message
    if fcbsupplied then fname = filename.filename else fname = filename
    message,'FITS file ' + fname + ' not modified',/INF
    return
   end 
