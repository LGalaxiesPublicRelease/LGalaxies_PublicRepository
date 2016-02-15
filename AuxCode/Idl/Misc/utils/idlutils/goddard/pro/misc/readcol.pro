pro readcol,name,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15, $
            v16,v17,v18,v19,v20,v21,v22,v23,v24,v25, COMMENT = comment, $
            FORMAT = fmt, DEBUG=debug, SILENT=silent, SKIPLINE = skipline, $
            NUMLINE = numline, DELIMITER = delimiter
;+
; NAME:
;       READCOL
; PURPOSE:
;       Read a free-format ASCII file with columns of data into IDL vectors 
; EXPLANATION:
;       Lines of data not meeting the specified format (e.g. comments) are 
;       ignored.  Columns may be separated by commas or spaces.
;
;       Use READFMT to read a fixed-format ASCII file.   Use RDFLOAT for
;       much faster I/O (but less flexibility).    Use FORPRINT to write 
;       columns of data (inverse of READCOL).    
;
; CALLING SEQUENCE:
;       READCOL, name, v1, [ v2, v3, v4, v5, ...  v25 , COMMENT=
;           DELIMITER= ,FORMAT = , /DEBUG ,  /SILENT , SKIPLINE = , NUMLINE = ]
;
; INPUTS:
;       NAME - Name of ASCII data file, scalar string.  
;
; OPTIONAL INPUT KEYWORDS:
;       FORMAT - scalar string containing a letter specifying an IDL type
;               for each column of data to be read.  Allowed letters are 
;               A - string data, B - byte, D - double precision, F- floating 
;               point, I - integer, L - longword, Z - longword hexadecimal, 
;               and X - skip a column.
;
;               Columns without a specified format are assumed to be floating 
;               point.  Examples of valid values of FMT are
;
;       'A,B,I'        ;First column to read as a character string, then 
;                       1 column of byte data, 1 column integer data
;       'L,L,L,L'       ;Four columns will be read as longword arrays.
;       ' '             ;All columns are floating point
;
;       If a FORMAT keyword string is not supplied, then all columns are 
;       assumed to be floating point.
;
;       /SILENT - Normally, READCOL will display each line that it skips over.
;               If SILENT is set and non-zero then these messages will be 
;               suppressed.
;       /DEBUG - If this keyword is non-zero, then additional information is
;                printed as READCOL attempts to read and interpret the file.
;       COMMENT - single character specifying comment signal.   Any line 
;                beginning with this character will be skipped.   Default is
;                no comment lines.
;       DELIMITER - single character specifying delimiter used to separate 
;                columns.   Default is either a comma or a blank.
;       SKIPLINE - Scalar specifying number of lines to skip at the top of file
;               before reading.   Default is to start at the first line.
;       NUMLINE - Scalar specifying number of lines in the file to read.  
;               Default is to read the entire file
;
; OUTPUTS:
;       V1,V2,V3,...V25 - IDL vectors to contain columns of data.
;               Up to 25 columns may be read.  The type of the output vectors
;               are as specified by FORMAT.
;
; EXAMPLES:
;       Each row in a file position.dat contains a star name and 6 columns
;       of data giving an RA and Dec in sexigesimal format.   Read into IDL 
;       variables.   (NOTE: The star names must not include the delimiter 
;       as a part of the name, no spaces or commas as default.)
;
;       IDL> FMT = 'A,I,I,F,I,I,F'
;       IDL> READCOL,'position.dat',F=FMT,name,hr,min,sec,deg,dmin,dsec  
;
;       The HR,MIN,DEG, and DMIN variables will be integer vectors.
;
;       Alternatively, all except the first column could be specified as
;       floating point.
;
;       IDL> READCOL,'position.dat',F='A',name,hr,min,sec,deg,dmin,dsec 
;
;       To read just the variables HR,MIN,SEC
;       IDL> READCOL,'position.dat',F='X,I,I,F',HR,MIN,SEC
;
; RESTRICTIONS:
;       This procedure is designed for generality and not for speed.
;       If a large ASCII file is to be read repeatedly, it may be worth
;       writing a specialized reader.
;
;       Columns to be read as strings must not contain the delimiter character
;       (i.e. commas or spaces by default).   Either change the default 
;       delimiter with the DELIMITER keyword, or use READFMT to read such files.
;
;       Numeric values are converted to specified format.  For example,
;       the value 0.13 read with an 'I' format will be converted to 0.
;
; PROCEDURES CALLED
;       GETTOK(), NUMLINES(), STRNUMBER()
;
; MINIMUM IDL VERSION:
;       V5.3 (Uses STRSPLIT)
; REVISION HISTORY:
;       Written         W. Landsman                 November, 1988
;       Modified             J. Bloch                   June, 1991
;       (Fixed problem with over allocation of logical units.)
;       Added SKIPLINE and NUMLINE keywords  W. Landsman    March 92
;       Read a maximum of 25 cols.  Joan Isensee, Hughes STX Corp., 15-SEP-93.
;       Call NUMLINES() function W. Landsman          Feb. 1996
;       Added DELIMITER keyword  W. Landsman          Nov. 1999
;       Fix indexing typos (i for k) that mysteriously appeared W. L. Mar. 2000
;       Hexadecimal support added.  MRG, RITSS, 15 March 2000.
;       Default is comma or space delimiters as advertised   W.L. July 2001
;       Faster algorithm, use STRSPLIT if V5.3 or later  W.L.  May 2002
;       Accept null strings separated by delimiter ,e.g. ',,,'
;       Use SCOPE_VARFETCH instead of EXECUTE() for >V6.1  W.L. Jun 2005
;       Added compile_opt idl2   W. L.  July 2005
;-
  On_error,2                           ;Return to caller
  compile_opt idl2

  if N_params() lt 2 then begin
     print,'Syntax - READCOL, name, v1, [ v2, v3,...v25, '
     print,'        FORMAT= ,/SILENT  ,SKIPLINE =, NUMLINE = , /DEBUG]'
     return
  endif

  no_exec = !VERSION.RELEASE GE '6.1'
; Get number of lines in file

   nlines = NUMLINES( name )
   if nlines LT 0 then return

   if keyword_set(DEBUG) then $
      message,'File ' + name+' contains ' + strtrim(nlines,2) + ' lines',/INF

   if not keyword_set( SKIPLINE ) then skipline = 0
   nlines = nlines - skipline
   if keyword_set( NUMLINE) then nlines = numline < nlines

  ncol = N_params() - 1           ;Number of columns of data expected
  vv = 'v' + strtrim( indgen(ncol)+1, 2)
  nskip = 0

  if N_elements(fmt) GT 0 then begin    ;FORMAT string supplied?

    if size(fmt,/tname) NE 'STRING' then $
       message,'ERROR - Supplied FORMAT keyword must be a scalar string'
;   Remove blanks from format string
    frmt = strupcase(strcompress(fmt,/REMOVE))   
    remchar, frmt, '('                  ;Remove parenthesis from format
    remchar, frmt, ')'           

;   Determine number of columns to skip ('X' format)
    pos = strpos(frmt, 'X', 0)

    while pos NE -1 do begin
        pos = strpos( frmt, 'X', pos+1)
        nskip = nskip + 1
    endwhile

  endif else begin                     ;Read everything as floating point

    frmt = 'F'
    if ncol GT 1 then for i = 1,ncol-1 do frmt = frmt + ',F'
    if not keyword_set( SILENT ) then message, $
      'Format keyword not supplied - All columns assumed floating point',/INF

  endelse

  nfmt = ncol + nskip
  idltype = intarr(nfmt)

; Create output arrays according to specified formats

   k = 0L                                     ;Loop over output columns
   hex = bytarr(nfmt)
   for i = 0L, nfmt-1 do begin

       fmt1 = gettok( frmt, ',' )
       if fmt1 EQ '' then fmt1 = 'F'         ;Default is F format
       case strmid(fmt1,0,1) of 
          'A':  idltype[i] = 7          
          'D':  idltype[i] = 5
          'F':  idltype[i] = 4
          'I':  idltype[i] = 2
          'B':  idltype[i] = 1
          'L':  idltype[i] = 3
          'Z':  begin 
                idltype[i] = 3               ;Hexadecimal
                hex[i] = 1b
                end
          'X':  idltype[i] = 0               ;IDL type of 0 ==> to skip column
         ELSE:  message,'Illegal format ' + fmt1 + ' in field ' + strtrim(i,2)
      endcase

; Define output arrays

      if idltype[i] GT 0 then begin
          if no_exec then (SCOPE_VARFETCH(vv[k], LEVEL=0))= $
	        make_array(nlines,TYPE = idltype[i]) else $
          tst = execute(vv[k] + '= make_array(nlines,TYPE = idltype[i] )' ) 
           k = k+1
      endif

   endfor
   goodcol = where(idltype)
   idltype = idltype[goodcol]
   check_numeric = (idltype NE 7)
   openr, lun, name, /get_lun
   ngood = 0L

   temp = ' '
   if !VERSION.RELEASE GE '5.6' then skip_lun,lun,skipline, /lines else $
   if skipline GT 0 then $
       for i = 0, skipline-1 do readf, lun, temp        ;Skip any lines

   if not keyword_set(delimiter) then delimiter = ' ,'
 
   for j = 0L, nlines-1 do begin

      readf, lun, temp
      if strlen(temp) LT ncol then begin    ;Need at least 1 chr per output line
          ngood = ngood-1
          if not keyword_set(SILENT) then $
                       message,'Skipping Line ' + strtrim(skipline+j+1,2),/INF
          goto, BADLINE 
       endif
    k = 0
    temp = strtrim(temp,1)                  ;Remove leading spaces
    if keyword_set(comment) then if strmid(temp,0,1) EQ comment then begin
          ngood = ngood-1
          if keyword_set(DEBUG) then $
                 message,'Skipping Comment Line ' + strtrim(skipline+j+1,2),/INF
          goto, BADLINE 
       endif

    var = strsplit(strcompress(temp),delimiter,/extract, /preserve_null) 
    if N_elements(var) LT nfmt then begin 
                 if not keyword_set(SILENT) then $ 
                      message,'Skipping Line ' + strtrim(skipline+j+1,2),/INF 
                 ngood = ngood-1            
                 goto, BADLINE         ;Enough columns?
    endif
    var = var[goodcol]

    for i = 0L,ncol-1 do begin
 
           if check_numeric[i] then begin    ;Check for valid numeric data
             tst = strnumber(var[i],val,hex=hex[i])          ;Valid number?
             if tst EQ 0 then begin            ;If not, skip this line
                 if not keyword_set(SILENT) then $ 
                      message,'Skipping Line ' + strtrim(skipline+j+1,2),/INF 
                 ngood = ngood-1
                 goto, BADLINE 
             endif
          if no_exec then $
	      (SCOPE_VARFETCH(vv[k], LEVEL=0))[ngood] = val else $
	       tst = execute(vv[k] + '[ngood] = val')

         endif else $
         if no_exec then $
	 (SCOPE_VARFETCH(vv[k], LEVEL=0))[ngood] = var[i] else $
           tst = execute(vv[k] + '[ngood] = var[i]')

      k = k+1

  endfor

BADLINE:  ngood = ngood+1

   endfor

  free_lun,lun
  if ngood EQ 0 then begin 
     message,'ERROR - No valid lines found for specified format',/INFORM
     return
  endif

  if not keyword_set(SILENT) then $
        message,strtrim(ngood,2) + ' valid lines read', /INFORM  

; Compress arrays to match actual number of valid lines
  if no_exec then begin
  for i=0,ncol-1 do $
       (SCOPE_VARFETCH(vv[i], LEVEL=0)) = $
            (SCOPE_VARFETCH(vv[i], LEVEL=0))[0:ngood-1]
 endif else begin
  for i = 0,ncol-1 do $
      tst = execute(vv[i] + '='+ vv[i]+ '[0:ngood-1]')
 endelse 

  return
  end
