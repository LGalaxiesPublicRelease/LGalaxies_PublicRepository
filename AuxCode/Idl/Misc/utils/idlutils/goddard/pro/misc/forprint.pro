pro forprint, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, $
      v15,v16,v17,v18,TEXTOUT = textout, FORMAT = format, SILENT = SILENT, $ 
      STARTLINE = startline, NUMLINE = numline, COMMENT = comment, $
      NoCOMMENT=Nocomment
;+
; NAME:
;       FORPRINT
; PURPOSE:
;       Print a set of vectors by looping over each index value.
;
; EXPLANATION:
;       If W and F are equal length vectors, then the statement
;               IDL> forprint, w, f   
;       is equivalent to 
;               IDL> for i = 0L, N_elements(w)-1 do print,w[i],f[i]    
;
; CALLING SEQUENCE:
;       forprint, v1,[ v2, v3, v4,....v18, FORMAT = , TEXTOUT = ,STARTLINE =,
;                                          NUMLINE =, /SILENT, COMMENT= ] 
;
; INPUTS:
;       V1,V2,...V18 - Arbitary IDL vectors.  If the vectors are not of
;               equal length then the number of rows printed will be equal
;               to the length of the smallest vector.   Up to 18 vectors
;               can be supplied.
;
; OPTIONAL KEYWORD INPUTS:
;
;       TEXTOUT - Controls print output device, defaults to !TEXTOUT
;
;               textout=1       TERMINAL using /more option if available
;               textout=2       TERMINAL without /more option
;               textout=3       file 'forprint.prt'
;               textout=4       file 'laser.tmp' 
;               textout=5      user must open file
;               textout =      filename (default extension of .prt)
;               textout=7       Append to <program>.prt file if it exists
;
;       COMMENT - String to write as the first line of output file if 
;                TEXTOUT > 2.    By default, FORPRINT will write a time stamp
;                on the first line.   Use /NOCOMMENT if you don't want FORPRINT
;                to write anything in the output file.
;       FORMAT - Scalar format string as in the PRINT procedure.  The use
;               of outer parenthesis is optional.   Ex. - format="(F10.3,I7)"
;               This program will automatically remove a leading "$" from
;               incoming format statements. Ex. - "$(I4)" would become "(I4)".
;               If omitted, then IDL default formats are used.
;       /NOCOMMENT  - Set this keyword if you don't want any comment line
;               line written as the first line in a harcopy output file.
;       /SILENT - Normally, with a hardcopy output (TEXTOUT > 2), FORPRINT will
;                print an informational message.    If the SILENT keyword
;               is set and non-zero, then this message is suppressed.
;       STARTLINE - Integer scalar specifying the first line in the arrays
;               to print.   Default is STARTLINE = 1, i.e. start at the
;               beginning of the arrays.
; OUTPUTS:
;       None
; SYSTEM VARIABLES:
;       If keyword TEXTOUT is not used, the default is the nonstandard 
;       keyword !TEXTOUT.    If you want to use FORPRINT to write more than 
;       once to the same file, or use a different file name then set 
;       TEXTOUT=5, and open and close then file yourself (see documentation 
;       of TEXTOPEN for more info).
;       
;       One way to add the non-standard system variables !TEXTOUT and !TEXTUNIT
;       is to use the procedure ASTROLIB
; EXAMPLE:
;       Suppose W,F, and E are the wavelength, flux, and epsilon vectors for
;       a spectrum.   Print these values to a file 'output.dat' in a nice 
;       format.
;
;       IDL> fmt = '(F10.3,1PE12.2,I7)'
;       IDL> forprint, F = fmt, w, f, e, TEXT = 'output.dat'
;
; PROCEDURES CALLED:
;       TEXTOPEN, TEXTCLOSE
; REVISION HISTORY:
;       Written    W. Landsman             April, 1989
;       Keywords textout and format added, J. Isensee, July, 1990
;       Made use of parenthesis in FORMAT optional  W. Landsman  May 1992
;       Added STARTLINE keyword W. Landsman    November 1992
;       Set up so can handle 18 input vectors. J. Isensee, HSTX Corp. July 1993
;       Handle string value of TEXTOUT   W. Landsman, HSTX September 1993
;       Added NUMLINE keyword            W. Landsman, HSTX February 1996
;       Added SILENT keyword             W. Landsman, RSTX, April 1998
;       Converted to IDL V5.0            W. Landsman, RSTX, April, 1998
;       Much faster printing to a file   W. Landsman, RITSS, August, 2001
;       Use SIZE(/TNAME) instead of DATATYPE() W. Landsman SSAI October 2001
;       Fix skipping of first line bug introduced Aug 2001  W. Landsman Nov2001
;       Added /NOCOMMENT keyword, the SILENT keyword now controls only 
;       the display of informational messages.  W. Landsman June 2002
;       Skip PRINTF if IDL in demo mode  W. Landsman  October 2004
;-            
  On_error,2                               ;Return to caller
  compile_opt idl2

  npar = N_params()
  if npar EQ 0 then begin
      print,'Syntax - FORPRINT, v1, [ v2, v3,...v18, FORMAT =, /SILENT, '
      print,'      /NoCOMMENT, COMMENT =, STARTLINE = , NUMLINE =, TEXTOUT =]'
      return
  endif

  if not keyword_set( STARTLINE ) then startline = 1l else $
         startline = startline > 1l 

  fmt="F"                 ;format flag
  istart = 2
  npts = N_elements(v1)

  if ( npts EQ 0 ) then message,'ERROR - Parameter 1 is not defined'

;  Remove "$" sign from format string and append parentheses if not 
;  already present

  if N_elements( format ) EQ 1 then begin

     fmt = "T"                                 ;format present
     frmt = format            
     if strmid(frmt,0,1) eq '$' then $
          frmt = strmid(frmt,1,strlen(frmt)-1) ;rem. '$' from format if present

     if strmid(frmt,0,1) NE '(' then frmt = '(' + frmt
     if strmid( frmt,strlen(frmt)-1,1) NE ')' then frmt = frmt + ')'

  endif

  if npar GT 1 then begin         ;Get number of elements in smallest array

      for i = istart, npar do begin 
          tst = execute('npts = min([ npts, N_elements(v'+strtrim(i,2)+')])')
          if npts EQ 0 then $
              message,'ERROR - Parameter ' + strtrim(i,2) + ' is not defined'
      endfor

  endif

  if keyword_set(NUMLINE) then npts = (startline + numline) < npts

  str = 'v' + strtrim( istart-1,2) + '[i]'
  if npar GT 1 then $
       for i = istart, npar do str = str + ',v' + strtrim(i,2) + '[i]'

; Use default output dev.
   demo = lmgr(/demo)
   if not demo then begin 

   if not keyword_set( TEXTOUT ) then textout = !TEXTOUT 
   if size( textout,/TNAME) EQ 'STRING' then text_out = 6  $      ;make numeric
                                  else text_out = textout

   textopen,'FORPRINT',TEXTOUT=textout,SILENT=silent
   if ( text_out GT 2 ) and (not keyword_set(NOCOMMENT)) then begin
       if N_elements(comment) GT 0 then $
        printf,!TEXTUNIT,comment else $
        printf,!TEXTUNIT,'FORPRINT: ',systime()
  endif 

; Determine whether we are printing to a terminal with /MORE capabilities.  If
; we are then we need to test for the user response at each pause.   Otherwise,
; we can print using a single EXECUTE call.

   more_term = (text_out EQ 1)          
   if more_term then begin
        stdout = fstat(-1)
        more_term =  stdout.isatty and (not stdout.isagui)
   endif
   endif
 
   if fmt EQ "F" then begin            ;Use default formats

   if demo then begin
         test =  execute('for i=startline-1,npts-1 do print,' + str)
        
   endif else if more_term then begin      
      for i = startline-1, npts-1 do begin 

          test = execute('printf,!TEXTUNIT,' + str) 
          if ( text_out EQ 1 ) then $                  ;Did user press 'Q' key?
               if !ERR EQ 1 then goto,DONE

      endfor
   endif else test = $
          execute('for i=startline-1,npts-1 do printf,!TEXTUNIT,' + str)

   endif else begin                    ;User specified format

   if demo then begin
         test =  execute('for i=startline-1,npts-1 do print,FORMAT=frmt,' + str)
 
   endif else  if more_term then begin

      for i = startline-1, npts-1 do begin 

         test = execute( 'printf, !TEXTUNIT,  FORMAT=frmt,' + str ) 
         if  ( text_out EQ 1 ) then  $
               if !ERR EQ 1 then goto,DONE

      endfor

    endif else test = $
        execute('for i=startline-1,npts-1 do printf,!TEXTUNIT,FORMAT=frmt,'+str)
        

  endelse

DONE: 
  textclose, TEXTOUT = textout          ;Close unit opened by TEXTOPEN

  return
  end
