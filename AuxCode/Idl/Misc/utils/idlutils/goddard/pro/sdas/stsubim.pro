pro stsubim, im, hd, filename, x1, x2, y1, y2, step, SILENT = silent
;+
; NAME:
;	STSUBIM
; PURPOSE:
;	Open an STSDAS file and read a portion of the file into an array.  
; EXPLANATION:
;	An enhanced version of STRD.  Program will prompt for the file name
;	and subimage bounds.
;
; CALLING SEQUENCE:
;	STSUBIM, im, hdr, [ filename, x1, x2, y1, y2, step, /SILENT ]
;
; OPTIONAL INPUTS:
;	FILENAME -  Character string giving the name of the SDAS file
;		to be read.  If omitted, then program will prompt 
;		for the file name.  If an extension is given, then
;		it must terminate in a 'h'.  A default extension of '.hhh' 
;		is assumed, if one is not supplied.  In VMS, version numbers 
;		are ignored, and the most recent version is always used.
;	X1      -  Starting x value (def=0), integer scalar
;	X2      -  Ending x value (def=total record size), integer scalar
;	Y1      -  Starting y value (def = 0), integer scalar
;	Y2      -  Ending y value (def = total no. of records), integer scalar
;     STEP    -  The number of pixels between extracted pixels.  This allows
;                the user to read every Nth pixel from the image.  STEP=1
;                indicates full resolution, STEP=2 indicates every other pixel,
;                etc.
; OUTPUTS:
;	IM - array containing image data
;	HDR - string array containing STSDAS header
;
; OPTIONAL KEYWORD INPUT:
;	SILENT  -  If this keyword is present, the size of the image 
;		will not be printed.
;
; COMMON BLOCKS:
;	STCOMMN - Created by SXOPEN.  STSUBIM uses STCOMMN to check
;	for an open unit, and to get image dimensions.          
;
; SIDE EFFECTS:
;	STSDAS image array and header are read into IM and HD
;	IF FILENAME is not supplied, then the program will check that
;	the image and header variable do not already contain data.
;
; RESTRICTIONS:
;	For use only on data without Groups!!  
;	For use only on 2 dimensional data files.
;
; PROCEDURE:
;	Program checks that STSDAS file exists and that IDL variables do
;	not already contain data, before calling SXOPEN and STSUB to
;	read in data.  The header keywords NAXIS* and CRPIX*
;	are updated to account for the actual image size.
;
; PROCEDURES CALLED:
;       FDECOMP, ORDINAL(), SPEC_DIR(), STSUB, SXOPEN, SXADDPAR, SXADDHIST
;
; MODIFICATION HISTORY:
;	Written B. Pfarr, STX, September 1987 from STRD
;	Modifed to IDL Version 2, M. Greason, STX, May 1990.
;	Prints 1st, 2nd, etc., instead of 1th, 2th... R. S.Hill, STX, Aug 1991
;	CRPIX transformation corrected.  RSH, HSTX, 27-May-1992.
;	Use new astrometry structure   W. Landsman   Feb 1994
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_error,2
 common stcommn,result,fname

 npar = N_params()
 err = string(7B) + 'STSUBIM: ERROR - '
 warn = string(7B) + 'STSUBIM: WARNING - '

;  Make sure the first three parameters, if supplied, are valid.

 if npar EQ 0 then begin
   print,'Syntax - STSUBIM, im,hdr, [ filename, x1, x2, y1, y2, step, /SILENT ]'
   return                                    
 endif
 
 if npar EQ 1 then begin
   ans = ''  
   print, warn,'A name for the header array was not supplied.'
   read, 'Continue the procedure to read the image array only [YES]? ',ANS
   if (strmid(strupcase(ans),0,1) EQ 'N') then return 
 endif

 if npar LE 2 then if ( N_elements(im) NE 0 ) then begin ;Already contain data?
   ans = ''
   print, warn, 'Image array contains data that will be erased'
   read,'Continue the procedure to read the image array [YES]? ',ans  
   if ( strmid( strupcase(ans),0,1 ) EQ 'N' ) then return 
 endif
 s = size(filename)

 if npar GE 3 then if (s[s[0] + 1] NE 7) then begin
	message,'Third parameter must be character string'
        return
	endif else file_name = filename  

 if NPAR LT 3 then begin
NAMER: 
   file_name = ''  ;Get file name if not supplied
   read, 'Enter name of SDAS data file (no quotes): ',file_name  
   if file_name EQ '' then return  
 endif

;  Check to see if file is on disk

 fdecomp, file_name, disk, dir, name, ext, ver
 if ver NE '' then file_name= disk + dir + name + '.' + ext ;No Versions allowed
 if ext EQ '' then file_name = file_name + '.hhh' $     ;Use default extension?
 else if strupcase(strmid(ext,2,1)) NE 'H' then begin
       	print, err, "SDAS file-name extensions must end with 'H'"
        goto, NAMER
 endif

 find = findfile( file_name, COUNT = ncount )  
 if ncount EQ 0 then begin
    print,err,'Unable to find ' + SPEC_DIR( file_name ) 
    if npar GE 3 then return
    print, 'Please re-enter name of file, or [RETURN] to exit'
    goto, NAMER
endif

;  Find an open unit between 1 and 9

 unit = 1    
 if N_elements(RESULT) NE 0 then begin	;Have other units been opened?
;	     If all 9 units have been used, then close unit 9 and reuse
 while ( result[1,unit] NE 0) and (unit LT 9) do unit = unit + 1  
 endif

;  Use SXOPEN to read header and initialize common block STCOMM

 sxopen, unit, file_name, hd
 xsiz = result[10,unit]  & ysiz = result[11,unit]

;  get subimage bounds and step

 if (npar NE 4) AND (npar LT 8) then begin
   step = ''
   read,'Enter step size of extraction (def =   1): ', step
   if step EQ '' then step = 1 else step = fix(step)
   step = step > 1
 endif else begin
   if npar EQ 4 then step = fix(x1)
 endelse

 if ( npar LE 4 ) then begin  
   x1 = ''
   read,'Enter lower limit for x (def =    0): ', x1  
   if x1 EQ '' then x1 = 0 else x1 = fix(x1)
 endif else x1 = fix(x1)

 if ( npar LT 5 ) then begin  
    x2 = ''
    read, 'Enter upper limit for x (def = '+ strtrim(xsiz-1,2) +'): ',x2 
    if x2 eq '' then x2 = xsiz-1 else x2=fix(x2)
 endif else x2 = fix(x2)

 if ( npar LT 6 ) then begin 
   y1 = ''
   read,'Enter lower limit for y (def =    0): ', Y1
   if y1 EQ '' then y1 = 0 else y1 = fix(y1)
 endif else y1 = fix(y1)

 if ( npar LT 7 ) then begin  
   y2=''
   read,'Enter upper limit for y (def = '+ strtrim(ysiz-1,2) +'): ',y2  
   if y2 eq '' then y2 = ysiz-1 else y2 = fix(y2)
 
 endif else y2 = fix(y2)

;  read in subimage

 if not keyword_set( SILENT ) then begin 
	message,'Now reading '+ strtrim((x2-x1)/step+1,2) +' by '+ $
	strtrim((y2-y1)/step+1,2) +' subarray',/INF
	if step NE 1 then message, $
                 'Reading every ' + ordinal(step) + ' pixel',/INF
 endif

 im = stsub( unit, x1, x2, y1, y2, step)
 close,unit

; Now adjust header for subimage

 sxaddpar,hd,'NAXIS1',(x2-x1+1)/step
 sxaddpar,hd,'NAXIS2',(y2-y1+1)/step

; Add history records

 label = 'STSUBIM:' + strmid(systime(),4,7) + strmid(systime(),20,4) +':'
 histrec = [label +  $
 ' Original Image Size Was '+ strtrim(xsiz,2) + ' by ' + strtrim(ysiz,2), $  
   label+' Original X Range: '+ strtrim(x1,2) +  ' to '+strtrim(x2,2),    $
   label+' Original Y Range: '+ strtrim(y1,2) +  ' to '+strtrim(y2,2)    ]
 if step NE 1 then histrec = [histrec, $
         label+' Every '+strupcase(ordinal(step))+ ' Original Pixel Used' ]
 sxaddhist, histrec, hd

; Update astrometry info if it exists

 extast, hd, astr, noparams

 if noparams LT 0 then return    ;If no astrometry in header then we are done

 if ( strmid(astr.ctype[0],5,3) EQ 'GSS') then begin
	gsss_stdast, hd
	extast, hd, astr, noparams
 endif

 crpix = astr.crpix
 sxaddpar, hd, 'CRPIX1', (crpix[0]-x1-0.5)/step + 0.5
 sxaddpar, hd, 'CRPIX2', (crpix[1]-y1-0.5)/step + 0.5

 if step EQ 1 then return     ;Plate scale unmodified unless STEP is non-unity

 if noparams EQ 2 then begin
	cd = astr.cd
	sxaddpar, hd, 'CD1_1', cd[0,0]*step
	sxaddpar, hd, 'CD1_2', cd[0,1]*step
	sxaddpar, hd, 'CD2_1', cd[1,0]*step
	sxaddpar, hd, 'CD2_2', cd[1,1]*step
 endif

 if (noparams EQ 0) or (noparams EQ 1) then begin
        cdelt = astr.cdelt
        sxaddpar, hd, 'CDELT1', cdelt[0]*step
        sxaddpar, hd, 'CDELT2', cdelt[1]*step
 endif      

 return
 end
