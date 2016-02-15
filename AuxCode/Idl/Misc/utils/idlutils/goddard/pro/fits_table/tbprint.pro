pro tbprint,hdr_or_tbstr,tab,columns,rows,textout=textout,fmt=fmt
;+
; NAME:
;       TBPRINT
;  PURPOSE:
;       Procedure to print specified columns & rows of a FITS binary table
;
; CALLING SEQUENCE:
;       TBPRINT, h, tab, columns, [ rows, TEXTOUT =, FMT = ]
;               or
;       TBPRINT,tb_str, tab, columns, [ rows, TEXTOUT =, FMT =  ]
;
; INPUTS:
;       h - FITS header for table, string array
;                       or
;       tb_str - IDL structure extracted from FITS header by TBINFO, useful 
;           when TBPRINT is called many times with the same header
;       tab - table array 
;       columns - string giving column names, or vector giving
;               column numbers (beginning with 1).  If string 
;               supplied then column names should be separated by comma's.
;       rows - (optional) vector of row numbers to print.  If
;               not supplied or set to scalar, -1, then all rows
;               are printed.
;
; OUTPUTS:
;       None
; OPTIONAL INPUT KEYWORDS:
;       TEXTOUT - scalar number (0-7) or string (file name) determining
;               output device (see TEXTOPEN).  Default is TEXTOUT=1, output 
;               to the user's terminal    
;       FMT = Format string for print display.   If not supplied, then any 
;               formats in the TDISP keyword fields of the table will be
;               used, otherwise IDL default formats.   
;
; SYSTEM VARIABLES:
;       Uses nonstandard system variables !TEXTOUT and !TEXTOPEN
;       Set !TEXTOUT = 3 to direct output to a disk file.   The system
;       variable is overriden by the value of the keyword TEXTOUT
;
; EXAMPLES:
;       tab = readfits('test.fits',htab,/ext) ;Read first extension into vars
;       tbprint,h,tab,'STAR ID,RA,DEC'    ;print id,ra,dec for all stars
;       tbprint,h,tab,[2,3,4],indgen(100) ;print columns 2-4 for 
;                                          first 100 stars
;       tbprint,h,tab,text="stars.dat"    ;Convert entire FITS table to
;                                         ;an ASCII file named 'stars.dat'
;
; PROCEDURES USED:
;       GETTOK(), TEXTOPEN, TEXTCLOSE, TBINFO
;
; RESTRICTIONS: 
;       (1) Program does not check whether output length exceeds output
;               device capacity (e.g. 80 or 132).
;       (2) Column heading may be truncated to fit in space defined by
;               the FORMAT specified for the column
;       (3) Program does not check for null values
;       (4) Does not work with variable length columns
;
; MINIMUM IDL VERSION:
;       V5.3 (uses STRSPLIT)
; HISTORY:
;       version 1  D. Lindler Feb. 1987
;       Accept undefined values of rows,columns W. Landsman  August 1997
;       Use new structure returned by TBINFO    W. Landsman  August 1997
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Made formatting more robust    W. Landsman   March 2000
;       Use STRSPLIT to parse string column listing W. Landsman July 2002
;       Wasn't always printing last row   W. Landsman  Feb. 2003
;-
 On_error,2


 if N_params() LT 2 then begin
   print,'Syntax -  TBPRINT, h, tab, [ columns, rows, device, TEXTOUT= ,FMT= ]'
   return
 endif

; set default parameters

 if N_elements(columns) EQ 0 then columns = -1
 if N_elements(rows) EQ 0 then rows= -1
 if not keyword_set(textout) then textout = 1
 nbytes = [1,2,4,4,8,8,1,0,16]
 fmt_def = ['','I4','I8','I12','G13.6','G16.8','','A','','','','']

; make sure rows is a vector

 sz = size(tab)
 nrows = sz[2]
 r = long(rows)
 if r[0] eq -1 then r = lindgen(nrows)          ;default
 n = N_elements(r)
;

 case  size(hdr_or_tbstr,/type) of 
 7: tbinfo,hdr_or_tbstr,tb_str
 8: tb_str = hdr_or_tbstr
 else: message,'ERROR - Invalid FITS header or structure supplied' 
 endcase 
 
 tfields = N_elements(tb_str.ttype)

; if columns is a string, change it to string array

 if size(columns,/tname) eq 'STRING' then begin
        colnames = strsplit(columns,',',/extract) 
        numcol = N_elements(colnames) 
        colnum = intarr(numcol)
        field = strupcase(colnames)
        for i = 0,numcol-1 do begin 
        colnum[i] = where(tb_str.ttype EQ field[i],nfound) + 1
        if nfound EQ 0 then $ 
           message,'Field '+ field[i] + ' not found in header'
       end
   end else begin                       ;user supplied vector
        colnum = fix(columns)           ;make sure it is integer
        numcol = N_elements(colnum)     ;number of elements
        if colnum[0] eq -1 then begin 
              colnum = indgen(tfields) + 1 & numcol = tfields
        endif 
 end

 if not keyword_set(fmt) then form = tb_str.tdisp[colnum-1] else begin
        if N_elements(fmt) EQ 1 and (numcol GT 1) then begin
                temp = strupcase(strtrim(fmt,2))
                if strmid(temp,0,1) EQ '(' then $
                        temp = strmid(temp,1,strlen(temp)-2)
                        form = strarr(numcol)
                        ifmt = 0
                         while strtrim(temp,2) NE ''  do begin
                                tstform = gettok(temp,',')
                                ndup = 1
                                vtype = strmid(tstform,0,1)
                                if strnumber(vtype,val) then begin
                                        ndup = val
                                        tstform = strmid(tstform,1,100)
                                endif
                                if strpos(tstform,'X') LT 0 then begin
                                     form[ifmt:ifmt+ndup-1]=tstform
                                     ifmt = ifmt + ndup
                                endif
                        endwhile
        endif else form = fmt
 endelse

 default = where(form EQ '',Ndef)
 if Ndef GT 0 then form[default] = fmt_def[ tb_str.idltype[colnum[default]-1] ]
 form = '(' + form + ')'
 
 num = where(tb_str.idltype[colnum-1] NE 7, Nnumeric)
 if Nnumeric GT 0 then minnumval = min(tb_str.numval[colnum[num]-1]) $
 else minnumval = 1
 if (minnumval GT 1) then begin 
        if rows[0] NE -1 then nrow1 = N_elements(rows)-1 else begin
                rows = lindgen(minnumval)
                nrow1 = minnumval-1
        endelse
        
 endif

 textopen,'TBPRINT', TEXTOUT = textout

 varname = 'v' + strtrim(sindgen(numcol)+1,2)
 len = lonarr(numcol)
 varstr = varname + '[0]'
 for i = 0,numcol-1 do begin
        result = execute(varname[i] + '= tbget(tb_str,tab,colnum[i],r)' )
        result = execute('len[i] = strlen(string(' + varstr[i] + ',f=form[i]))')
 endfor

 field = strtrim(tb_str.ttype[colnum-1],2)
 fieldlen = strlen(field)
 for i=0,numcol-1 do begin
        if fieldlen[i] LT len[i]-1 then begin
                pad = string(replicate(32b,(len[i] - fieldlen[i])/2))
                field[i] = pad + field[i] + pad
        endif else field[i] = strmid(field[i],0,len[i])
 endfor

 printf,!TEXTUNIT,field

 varstr = varname + '[i]'
 vstring = varstr[0]
 if numcol GT 1 then for i=1,numcol-1 do vstring = vstring + ',' + varstr[i]
 format = form[0]
 if N_elements(form) GT 0 then for i=1,N_elements(form)-1 do $
        format = format + ',' + form[i]
 format = '(' + format + ')'


 if minnumval EQ 1 then $
 result = execute('for i=0,n-1 do printf,!TEXTUNIT,' +  $
                   vstring + ',f=format') else $
 result = execute('for i=rows[0],rows[nrow1] do printf,!TEXTUNIT,' +  $
                   vstring + ',f=fmt') 
 textclose, TEXTOUT = textout
 return
 end
