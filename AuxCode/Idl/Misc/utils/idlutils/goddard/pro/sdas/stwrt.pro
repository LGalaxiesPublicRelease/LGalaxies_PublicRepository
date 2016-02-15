pro stwrt, im, hd, name, SDAS=sdas, NOUPDATE = noupdate
;+
; NAME:
;	STWRT
; PURPOSE:
;	Write a STSDAS header and/or 2-D array to disk (without groups!) 
;
; CALLING SEQUENCE:
;	STWRT, hdr                     ;Write image header (.hhh file) only
;	STWRT, im                      ;Create header to match image array
;	STWRT, im, hdr,[ name, /SDAS, /NOUPDATE ]  
;
; INPUT PARAMETERS:
;	im - image array to be written to disk.  If no header array is
;		supplied, then a simple header appropiate to the image will be
;		constructed.
;
; OPTIONAL INPUT PARAMETER:
;	hdr - STSDAS header, string array.  
;	name - character string containing the name of output file
;		to which the image is written.  If omitted, then 
;		the program will prompt for the file name.  A file
;		name will have the default extension of '.HHH'
;
; OPTIONAL KEYWORD INPUTS:
;	NOUPDATE-- By default, STWRT will modify the input FITS header to create
;		a proper SDAS .hhh file.   This includes  ensuring that (1) a 
;               DATATYPE keyword exists, and (2) that BITPIX is a positive
;		value.   
;	SDAS-- The SDAS keyword can be specified to modify the input header to 
;		ensure SDAS compatibility.
;
; RESTRICTIONS:
;	(1) STWRT is only for 2 dimensional images.  For other arrays use
;
;	SXOPEN,1,NAME,HD,HISTORY,'W'   ;HISTORY need not be defined
;	SXWRITE,1,IM
;	CLOSE,1
;
;	(2) The type of data written is determined by the DATATYPE 
;		keyword in the header.  A DATATYPE keyword appropiate to
;		the image array type will be added if does not already exist
;
; SIDE EFFECTS:
;	A STSDAS header and/or image array is written to disk.    Unit 2 is
;	opened and closed.
;
; REVISION HISTORY:
;	Written W. Landsman, STI Corp. August, 1986
;	Returned old version to not modify header.  W. Landsman July 1991.
;	Included call to CHKDType--option to change DATATYPE and BITPIX if they; 
;	do not match.				J.D.Offenberg Dec 1991.
;	Added call to CHECK_FITS, NOUPDATE keyword, remove autochange keyword
;                                        W. Landsman May 1992
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
 npar = N_params() 

 case npar of 
 0: begin
   print,'Syntax: STWRT, im, hdr, [ filename]   ;Write image & header' 
   print,'    or: STWRT, im                   ;Write image with minimal header
   print,'    or: STWRT, hdr                  ;Write header only'
   return
   end

 1: begin
   info = size(im)
   if (info[2] EQ 7) then begin		;Header supplied-first parameter?
        name = ''
        read,'Enter name of output file (no quotes): ',name
	sxhwrite, name, im
        return
   endif else begin
   Ans = 'Y'
   print,'No header array supplied'
   print,'Do you want to create a minimal new header based on the
   read,'information in the supplied array [YES]: ', Ans
   if strlen(ans) EQ 0 then ans = 'Y'
   If strupcase(strmid(ans,0,1)) NE 'Y' then begin
        print,'Program returning: no files created'
        return
   endif
   endelse
   end

 else: if keyword_set(NOUPDATE) then check_FITS, im, hd  $
          else check_FITS, im, hd, /UPDATE
 endcase
 
 if npar LE 2 then begin	;Prompt for file name?
     name = ''
     read,'Enter name of output file (no quotes): ',name
 endif

 close,2

 if (npar gt 1) then if keyword_set(SDAS) then $
           sxmake, 2, name, im, 0, 0, hd else $   ;write image and header
           sxopen, 2, name, hd, history, 'W' $    ;write image and header
      else sxmake, 2, name, im, 0, 0                    ;write image only

 sxwrite, 2, im
 close, 2

 return
 end
