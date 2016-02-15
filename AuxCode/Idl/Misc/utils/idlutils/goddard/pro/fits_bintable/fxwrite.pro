	PRO FXWRITE, FILENAME, HEADER, DATA, NANVALUE=NANVALUE,		$
		NOUPDATE=NOUPDATE, ERRMSG=ERRMSG
;+
; NAME: 
;	FXWRITE
; Purpose     : 
;	Write a disk FITS file.
; Explanation : 
;	Creates a disk FITS file and writes a FITS primary header, and
;	optionally a primary data array.
; Use         : 
;	FXWRITE, FILENAME, HEADER [, DATA ]
; Inputs      : 
;	FILENAME = String containing the name of the file to be written.
;	HEADER	 = String array containing the header for the FITS file.
; Opt. Inputs : 
;	DATA	 = IDL data array to be written to the file.  If not passed,
;		   then it is assumed that extensions will be added to the
;		   file.
; Outputs     : 
;	None.
; Opt. Outputs: 
;	None.
; Keywords    : 
;	NANVALUE = Value signalling data dropout.  All points corresponding to
;		   this value are set to be IEEE NaN (not-a-number).  Ignored
;		   unless DATA is of type float, double-precision or complex.
;	NOUPDATE = If set, then the optional BSCALE and BZERO keywords in the
;		   HEADER array will not be changed.  The default is to reset
;		   these keywords to BSCALE=1, BZERO=0.
;	ERRMSG	 = If defined and passed, then any error messages will be
;		   returned to the user in this parameter rather than
;		   depending on the MESSAGE routine in IDL.  If no errors are
;		   encountered, then a null string is returned.  In order to
;		   use this feature, ERRMSG must be defined first, e.g.
;
;			ERRMSG = ''
;			FXWRITE, ERRMSG=ERRMSG, ...
;			IF ERRMSG NE '' THEN ...
;
; Calls       : 
;	CHECK_FITS, GET_DATE, HOST_TO_IEEE, FXADDPAR, FXPAR
; Common      : 
;	None.
; Restrictions: 
;	If DATA is passed, then HEADER must be consistent with it.  If no data
;	array is being written to the file, then HEADER must also be consistent
;	with that.  The routine FXHMAKE can be used to create a FITS header.
;
;	If found, then the optional keywords BSCALE and BZERO in the HEADER
;	array is changed so that BSCALE=1 and BZERO=0.  This is so that these
;	scaling parameters are not applied to the data a second time by another
;	routine.  Also, history records are added storing the original values
;	of these constants.  (Other values of BZERO are used for unsigned
;	integers.)
;
;	If the /NOUPDATE keyword is set, however, then the BSCALE and BZERO
;	keywords are not changed.  The user should then be aware that FITS
;	readers will apply these numbers to the data, even if the data is
;	already converted to floating point form.
;
;	Groups are not supported.
;
; Side effects: 
;	None.
; Category    : 
;	Data Handling, I/O, FITS, Generic.
; Prev. Hist. : 
;	W. Thompson, Jan 1992, from WRITEFITS by J. Woffard and W. Landsman.
;	Differences include:
;
;		* Made DATA array optional, and HEADER array mandatory.
;		* Changed order of HEADER and DATA parameters.
;		* No attempt made to fix HEADER array.
;
;	W. Thompson, May 1992, changed open statement to force 2880 byte fixed
;			       length records (VMS).  The software here does not
;			       depend on this file configuration, but other
;			       FITS readers might.
;	W. Thompson, Aug 1992, added code to reset BSCALE and BZERO records,
;			       and added the NOUPDATE keyword.
; Written     : 
;	William Thompson, GSFC, January 1992.
; Modified    : 
;	Version 1, William Thompson, GSFC, 12 April 1993.
;		Incorporated into CDS library.
;	Version 2, William Thompson, GSFC, 31 May 1994
;		Added ERRMSG keyword.
;	Version 3, William Thompson, GSFC, 23 June 1994
;		Modified so that ERRMSG is not touched if not defined.
;	Version 4, William Thompson, GSFC, 12 August 1999
;		Catch error if unable to open file.
;       Version 4.1 Wayne Landsman, GSFC, 02 May 2000
;               Remove !ERR in call to CHECK_FITS, Use ARG_PRESENT()
;       Version 5, William Thompson, GSFC, 22 September 2004
;               Recognize unsigned integer types
; Version     : 
;	Version 5, 22 Sep 2004
;-
;
	ON_ERROR, 2
;
;  Check the number of parameters.
;   
	IF N_PARAMS() LT 2 THEN BEGIN
	    MESSAGE = 'Syntax:  FXWRITE, FILENAME, HEADER  [, DATA ]'
	    GOTO, HANDLE_ERROR
	ENDIF
;
;  Check the header against the data being written to the file.  If the data
;  array is not passed, then NAXIS should be set to zero, and EXTEND should be
;  true.
;
	IF N_PARAMS() EQ 2 THEN BEGIN
	    IF (FXPAR(HEADER,'NAXIS') NE 0) THEN BEGIN
		MESSAGE = 'NAXIS should be zero for no primary data array'
		GOTO, HANDLE_ERROR
	    END ELSE IF (NOT FXPAR(HEADER,'EXTEND')) THEN BEGIN
		MESSAGE = 'EXTEND should be true for no primary data array'
		GOTO, HANDLE_ERROR
	    ENDIF
	END ELSE BEGIN
	    CHECK_FITS, DATA, HEADER, /FITS, ERRMSG = MESSAGE
	    IF MESSAGE NE '' THEN GOTO, HANDLE_ERROR
	ENDELSE
;
;  Set the BSCALE and BZERO keywords to their default values.
;
        SZ = SIZE(DATA)
        TYPE = SZ[SZ[0]+1]
        IF N_PARAMS() EQ 3 THEN NEWDATA = DATA
	IF NOT KEYWORD_SET(NOUPDATE) THEN BEGIN
	    BZERO  = FXPAR(HEADER,'BZERO')
	    BSCALE = FXPAR(HEADER,'BSCALE')
	    GET_DATE,DTE
	    IF (BSCALE NE 0) AND (BSCALE NE 1) THEN BEGIN
		FXADDPAR,HEADER,'BSCALE',1.
		FXADDPAR,HEADER,'HISTORY',DTE+' reset BSCALE, was '+ $
			STRTRIM(BSCALE,2)
            ENDIF
;
;  If an unsigned data type then redefine BZERO to allow all the data to be
;  stored in the file.
;
            BZERO0 = 0
            IF (TYPE EQ 12) AND (NOT KEYWORD_SET(NOUPDATE)) THEN BEGIN
                BZERO0 = '8000'X
                NEWDATA = FIX(TEMPORARY(NEWDATA) - BZERO)
            ENDIF
            IF (TYPE EQ 13) AND (NOT KEYWORD_SET(NOUPDATE)) THEN BEGIN
                BZERO0 = '80000000'X
                NEWDATA = LONG(TEMPORARY(NEWDATA) - BZERO)
            ENDIF
	    IF BZERO NE BZERO0 THEN BEGIN
		FXADDPAR,HEADER,'BZERO',BZERO0
		FXADDPAR,HEADER,'HISTORY',DTE+' reset BZERO, was '+ $
			STRTRIM(BZERO,2)
	    ENDIF
	ENDIF
;
;  Get the UNIT number, and open the file.
;
       	GET_LUN, UNIT      
       	OPENW, UNIT, FILENAME, 2880, /BLOCK, ERROR=ERR
	IF ERR NE 0 THEN BEGIN
	    MESSAGE = 'Error creating file '+FILENAME
	    GOTO, HANDLE_ERROR
	ENDIF
;
;  Determine if an END line occurs, and add one if necessary
;
	ENDLINE = WHERE( STRMID(HEADER,0,8) EQ 'END     ', NEND)
	ENDLINE = ENDLINE[0]
	IF NEND EQ 0 THEN BEGIN
	    MESSAGE, 'WARNING - An END statement has been appended ' + $
		'to the FITS header', /INFORMATIONAL
	    HEADER = [HEADER, 'END' + STRING(REPLICATE(32B,77))]
	    ENDLINE = N_ELEMENTS(HEADER) - 1 
	ENDIF
	NMAX = ENDLINE + 1		;Number of 80 byte records
	NHEAD = FIX((NMAX+35)/36)	;Number of 2880 byte records
;
;  Convert to byte and force into 80 character lines
;
	BHDR = REPLICATE(32B, 80, 36*NHEAD)
	FOR N = 0,ENDLINE DO BHDR[0,N] = BYTE( STRMID(HEADER[N],0,80) )
	WRITEU, UNIT, BHDR
;
;  If passed, then write the data array.
;
	IF N_PARAMS() EQ 3 THEN BEGIN
;
;  If necessary, then byte-swap the data before writing it out.  Also, replace
;  any values corresponding data dropout with IEEE NaN.
;
	    IF (N_ELEMENTS(NANVALUE) EQ 1) AND (TYPE GE 4) AND	$
		    (TYPE LE 6) THEN BEGIN
		W = WHERE(DATA EQ NANVALUE, COUNT)
		CASE TYPE OF
		    4:  NAN = FLOAT(  REPLICATE('FF'XB,4),0,1)
		    5:  NAN = DOUBLE( REPLICATE('FF'XB,8),0,1)
		    6:  NAN = COMPLEX(REPLICATE('FF'XB,8),0,1)
		    9:  NAN = DCOMPLEX(REPLICATE('FF'XB,16),0,1)
		ENDCASE
	    END ELSE COUNT = 0
;
	    HOST_TO_IEEE, NEWDATA
	    IF COUNT GT 0 THEN NEWDATA[W] = NAN
;
	    WRITEU,UNIT,NEWDATA
;
;  If necessary, then pad out to an integral multiple of 2880 bytes.
;
	    BITPIX = FXPAR( HEADER, 'BITPIX' )
	    NBYTES = N_ELEMENTS(DATA) * (ABS(BITPIX) / 8 )
	    NPAD = NBYTES MOD 2880
	    IF NPAD NE 0 THEN BEGIN
		NPAD = 2880 - NPAD
		WRITEU,UNIT,BYTARR(NPAD)
	    ENDIF
	ENDIF
;
;  Close the file and return.
;
	FREE_LUN, UNIT
	IF ARG_PRESENT(ERRMSG)  THEN ERRMSG = ''
	RETURN
;
HANDLE_ERROR:
	IF N_ELEMENTS(UNIT) EQ 1 THEN FREE_LUN, UNIT
	IF ARG_PRESENT(ERRMSG) THEN ERRMSG = 'FXWRITE: ' + MESSAGE	$
		ELSE MESSAGE, MESSAGE
;
	END
