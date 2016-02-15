pro ftab_ext,file_or_fcb,columns,v1,v2,v3,v4,v5,v6,v7,v8,v9,ROWS=rows, $
        EXTEN_NO = exten_no
;+
; NAME:
;       FTAB_EXT
; PURPOSE:
;       Routine to extract columns from a FITS (binary or ASCII) table. 
;
; CALLING SEQUENCE:
;       FTAB_EXT, name_or_fcb, columns, v1, [v2,..,v9, ROWS=, EXTEN_NO= ]
; INPUTS:
;       name_or_fcb - either a scalar string giving the name of a FITS file 
;               containing a (binary or ASCII) table, or an IDL structure 
;               containing as file control block (FCB) returned by FITS_OPEN 
;               If FTAB_EXT is to be called repeatedly on the same file, then
;               it is quicker to first open the file with FITS_OPEN, and then
;               pass the FCB structure to FTAB_EXT
;       columns - table columns to extract.  Can be either 
;               (1) String with names separated by commas
;               (2) Scalar or vector of column numbers
;
; OUTPUTS:
;       v1,...,v9 - values for the columns.   Up to 9 columns can be extracted
;
; OPTIONAL INPUT KEYWORDS:
;       ROWS -  scalar or vector giving row number(s) to extract
;               Row numbers start at 0.  If not supplied or set to
;               -1 then values for all rows are returned
;       EXTEN_NO - Extension number to process.   If not set, then data is
;               extracted from the first extension in the file (EXTEN_NO=1)
;
; EXAMPLES:
;       Read wavelength and flux vectors from the first extension of a 
;       FITS file, 'spec.fit'.   Using FTAB_HELP,'spec.fit' we find that this
;       information is in columns named 'WAVELENGTH' and 'FLUX' (in columns 1
;       and 2).   To read the data
;
;       IDL> ftab_ext,'spec.fit','wavelength,flux',w,f
;               or
;       IDL> ftab_ext,'spec.fit',[1,2],w,f
;       
; PROCEDURES CALLED:
;       FITS_READ, FITS_CLOSE, FTINFO, FTGET(), TBINFO, TBGET()
; MINIMUM IDL VERSION:
;       V5.3 (uses STRSPLIT)
; HISTORY:
;       version 1        W.   Landsman         August 1997
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Improve speed processing binary tables  W. Landsman   March 2000
;       Use new FTINFO calling sequence  W. Landsman   May 2000  
;       Don't call fits_close if fcb supplied W. Landsman May 2001 
;       Use STRSPLIT to parse column string  W. Landsman July 2002 
;       Cleanup pointers in TBINFO structure  W. Landsman November 2003
;-
;---------------------------------------------------------------------
 if N_params() LT 3 then begin
        print,'Syntax - FTAB_EXT, name, columns, v1, [v2,...,v9, ROWS=, EXTEN=]'
        return
 endif
 N_ext = N_params() - 2
 strng = size(columns,/TNAME) EQ 'STRING'    ;Is columns a string?

 if not keyword_set(exten_no) then exten_no = 1
 sz = size(file_or_fcb)
 if sz[sz[0]+1] NE 8 then fits_open,file_or_fcb,fcb else fcb=file_or_fcb
 if fcb.nextend EQ 0 then $
        message,'ERROR - FITS file contains no table extensions'
 if fcb.nextend LT exten_no then $
        message,'ERROR - FITS file contains only ' + strtrim(fcb.nextend,2) + $
                ' extensions'

 if N_elements(rows) NE 0 then begin
        minrow = min(rows, max = maxrow)
        naxis1 = fcb.axis[0,exten_no]
        first = naxis1*minrow
        last = naxis1*(maxrow+1)-1
        xrow = rows - minrow
        fits_read,fcb,tab,htab,exten_no=exten_no,first=first,last=last,/no_pdu
        tab = reform(tab,naxis1,maxrow-minrow+1,/overwrite)
 endif else begin
        fits_read, fcb, tab, htab, exten_no=exten_no,/no_pdu 
        xrow = -1
 endelse
 if sz[sz[0]+1] NE 8 then fits_close,fcb else $
         file_or_fcb.last_extension = exten_no
 ext_type = fcb.xtension[exten_no]

 case ext_type of
 'A3DTABLE': binary = 1b
 'BINTABLE': binary = 1b
 'TABLE': binary = 0b
 else: message,'ERROR - Extension type of ' + $
                ext_type + 'is not a FITS table format'
 endcase

 if strng then colnames= strsplit(columns,',',/EXTRACT) else $
               colnames = columns
 if binary then tbinfo,htab,tb_str else ftinfo,htab,ft_str

 for i = 0, N_ext-1 do begin
         if binary then $
                v = TBGET( tb_str,tab,colnames[i],xrow,nulls) $
        else $
                v = FTGET( ft_str,tab,colnames[i],xrow,nulls)
        command = 'v'+strtrim(i+1,2)+'=v'
        istat = execute(command)
        endfor
 if binary then begin
        ptr_free, tb_str.tscal
        ptr_free, tb_str.tzero
 endif
 return
 end 


