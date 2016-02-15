pro strd, im, hd, filename, GROUP = group, PAR = PAR 
;+
; NAME:
;	STRD
; PURPOSE:
;	Open an STSDAS file and read into an image array and header. 
; EXPLANATION:
;	Combines the functions of SXREAD and SXOPEN. 
;
; CALLING SEQUENCE:
;	STRD, im, hdr, [ filename, GROUP = , PAR = ]
;
; OPTIONAL INPUT:
;	FILENAME -  Character string giving the name of the SDAS file
;		to be read.  If omitted, then program will prompt 
;		for the file name.  If an extension is given, then
;		it must terminate in a 'h'.
;		A default extension of '.hhh' is assumed, if one is
;		not supplied.  Under VMS, the version numbers are ignored, 
;		and the most recent version is always used.
;
; OUTPUTS:
;	IM - array containing image data
;	HDR - string array containing header
;
; OPTIONAL INPUT KEYWORD PARAMETER:
;	GROUP - scalar integer specifying group number to read.  Default is 0.
;
; OPTIONAL OUTPUT KEYWORD PARAMETER:
;	PAR - Parameter block (byte array) read from group formatted data
; COMMON BLOCKS:
;	STCOMMN - Created by SXOPEN.  STRD uses STCOMMN to check
;		for an open unit, and to get image dimensions.          
;
; SIDE EFFECTS:
;	STSDAS image array and header are read into IM and HDR
;	IF FILENAME is not supplied, then the program will check that
;	the image and header variable do not already contain data.
;
; SYSTEM VARIABLES:
;	Set !QUIET = 1 to suppress informational messages.
;
; METHOD:
;	Program checks that specified STSDAS file exists before calling 
;	SXOPEN and SXREAD to read in data.
;
; PROCEDURES CALLED:
;	FDECOMP, PICKFILE(), SPEC_DIR(), SXOPEN, SXREAD()
; MODIFICATION HISTORY:
;	Written W. Landsman, STI Corporation August 1986
;	Optional parameter "FILENAME" added November 1986
;	Correctly print header size when more than 2 dimensions February 1996
;	Add GROUP, PAR keywords, call PICKFILE   W. Landsman   March 1996
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_error,2

 common stcommn,result,fname
 npar = N_params()
 
 if npar EQ 0 then begin
   print,'Syntax - STRD, im, hd, [filename, GROUP =, PAR = ]
   return
 endif

 if not keyword_set(GROUP) then group = 0

  if ( N_elements( filename ) EQ 1) then  begin 
		file_name = filename
		fdecomp, file_name, disk, dir, name, ext
	      	if name NE '' then goto, FINDER
 endif 

NAMER:

 if ( !D.FLAGS AND 65536) NE 0 then begin       ;Have widgets?
     if N_elements(dir) EQ 1 then $
          file_name = pickfile(filter = '*.hhh',/NOCON,/READ, $
		path = disk + dir) else $
          file_name = PICKFILE(filter = '*.hhh',/NOCON,/Read)
 endif else begin   
   file_name = ''          ;Get file name if not supplied
   read,'Enter name of SDAS data file (no quotes): ',file_name  
 endelse
 if ( file_name  EQ '') then return 

FINDER: 
 fdecomp, file_name, disk, dir, name, ext, ver   

 if ver NE '' then $
    file_name = disk + dir + name + '.' + ext          ;No Versions allowed

 if ext EQ '' then file_name = file_name + '.hhh' $     ;Use default extension?
 else if strupcase( strmid(ext,2,1) ) NE 'H' then begin     
      	message, "ERROR - SDAS file-name extensions must end with 'h'",/CON
        goto, NAMER 
 endif

 find = findfile( file_name, COUNT = I )    ;Does file exist?
 if I LT 1 then begin
    message,'ERROR - Unable to find ' + SPEC_DIR( file_name ), /CON
    if npar EQ 3 then return   
    print, 'Please re-enter name of file, or [RETURN] to exit'
    GOTO, namer  
 endif 

 for i = 1, 9 do begin               ;Find an open unit between 1 and 9
   test = fstat(i)
   if test.open eq 0 then begin
      unit = i
      goto, OPENER      
   endif
 endfor

OPENER: 

  sxopen, unit, file_name, hd   
  naxis = result[3,unit]
  dimen = strtrim(result[10:10+naxis-1,unit],2)
  st = dimen[0]
  if naxis GT 1 then for i=1,naxis-1 do st = st + ' by ' + dimen[i] $
	else st = st + ' element'

  message, /INF, 'Now reading '+ st + ' array'

 im = sxread(unit,group,par) 

 close, unit    
 return
 end
