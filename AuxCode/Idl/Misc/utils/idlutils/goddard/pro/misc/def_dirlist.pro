	PRO DEF_DIRLIST, EVAR, VALUE
;+
; NAME:	
;       DEF_DIRLIST
;
; PURPOSE:
;       Define directory list using setenv or setlog
;
; EXPLANATION:	
;       Environment variables which point to a list of directories can
;       end up to be very long.  In VMS this can be a problem, because logical 
;       names cannot be longer than 256 characters.  However, it is possible to
;       get around this in VMS by assigning multiple values to a single logical
;       name--a facility that does not exist in Unix.
;
;       This routine will define the environment variable as either a delimited
;       string, or as a series of values, whichever is most appropriate.
;
; CALLING SEQUENCE:	
;       DEF_DIRLIST, EVAR, VALUE
; INPUTD:	
;       EVAR = The name of the environment variable to define.
;       VALUE = The value to give to EVAR.  This can be either a single, 
;               delimited string, or it can be an array of directory names.
;               The routine will choose for itself how to use this to define the 
;               environment variable.
;
; EXAMPLES:	
;       DIRS = FIND_ALL_DIR('+/data/fits')
;       DEF_DIRLIST, 'FITS_DATA', DIRS
;
; PROCEDURE CALLS:
;       SETENV, STR_SEP()
;	Note: The intrinsic SETENV command is available under Unix & Windows
;	only.   However, it is available as a Library procedure for VMS.
;
; REVISION HISTORY:	
;	Version 1, 06-Aug-1996, William Thompson, GSFC
;       Converted to IDL V5.0   June 1998    W. Landsman
;       Use STRSPLIT instead of STR_SEP if V5.3 or later W.L.  July 2002
;-
;
	ON_ERROR, 2
;
	IF N_PARAMS() NE 2 THEN MESSAGE, 'Syntax:  DEF_DIRLIST, EVAR, VALUE'
;
;  Form a delimited string from the input.
;
	DIR = VALUE[0]
	IF !VERSION.OS EQ 'vms' THEN SEP = ',' ELSE SEP = ':'
	FOR I = 1,N_ELEMENTS(VALUE)-1 DO DIR = DIR + SEP + VALUE[I]
;
;  If VMS, and the string length is greater than 256, then use SETLOG to define
;  the value.
;
	IF !VERSION.OS_FAMILY EQ 'vms' AND STRLEN(DIR) GE 256 THEN BEGIN
                   DIR = STRSPLIT(DIR,SEP,/EXTRACT) 
		SETLOG, EVAR, DIR
;
;  Otherwise, use SETENV
;
	END ELSE SETENV, STRTRIM(EVAR,2) + '=' + STRTRIM(DIR,2)
;
	RETURN
	END
