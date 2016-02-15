PRO MID_UP_IMAGE,IMAGE,DATA,NAXIS,NPIX
;-------------------------------------------------------------------------------
;+
; NAME:
;	MID_UP_IMAGE
;
; PURPOSE:
;	Get a pixel matrix and some support information from a MIDAS file.
; EXPLANATION:
;	Allows updating of data, using DATA parameter.
;
; CALLING SEQUENCE:
;	MID_UP_IMAGE,IMAGE,DATA,NAXIS,NPIX
;
; INPUTS:
;	IMAGE = Filename or Logical Unit Number.
;	* If a filename is given, the file will be opened and closed using a
;		local LUN.  The filename is that of the MIDAS image, without
;		extension (.BDF will is assumed) or version number (latest 
;		version is assumed).
;	* If a LUN is given, the file associated with that LUN will be used.
;
; OUTPUTS:
;	NAXIS = Number of dimensions in MIDAS image. I*4 values.
;	NPIX  = Array containing the dimensions of the data to be written into
;		the MIDAS image. Must be compatible with (i.e. smaller than or
;		equal to the corresponding dimensions of) the latter. 
;		I*4 values.
;	DATA  = Array to be written into the MIDAS image. Dimensions are 
;		defined by NAXIS and NPIX. R*4 values.
;
; ALGORITHM:
;	0) Check inputs and set error handling.
;	1) Open file for access using the access method indicated by the type of
;		the input parameter IMAGE.
;	2) Get the descriptors of the IMAGE.
;		a) NAXIS set the output parameter NAXIS to this value.
;		b) NPIX  set the output parameter NPIX to this value.
;	3) Locate the pixel data start block
;	4) Load pixel data into output parameter DATA
;	5) Check consistency of dimensions of data array to be written, and
;		image dimensions; then write data array into image.
;	6) Terminate file access as is proper for the type of parameter IMAGE.
;
; RESTRICTIONS:
;	1) There must be only one FCB and it must be in VBN 1.
;	2) The LDBs must begin in VBN 2.
;	3) All descriptors must be in the first LDB.
;	4) Currently only works for real data, does not check to see if this is
;		true or not.
;	5) Midas extensions (.bdf, .tbl) assumed lower case.
;
; AUTHORS:
;	FM   - F. Murtagh, ST-ECF
;	SAV  - Stephen A. Voels, USM/DAN
;
; MODIFICATION HISTORY:
;	MAY 1989  FM   Initial programming.
;	FEB 1991  FM   Conversion to v.2 IDL, Unix.
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-------------------------------------------------------------------------------
; STEP 0:
 
; Set defaults for values to be derived from the MIDAS image.
 
NAXIS_IM =  0L
NPIX_IM  = [0L]
DAT  = [0.0]
 
; Verify the number of parameters and the types of the inputs.
 
N = N_PARAMS(0)
IF (N NE 4) THEN BEGIN
   PRINT,'CALLING SEQUENCE: MID_UP_IMAGE,IMAGE,DATA,NAXIS,NPIX'
   RETURN
   ENDIF
 
INFO_IMAGE = SIZE(IMAGE)
IF (INFO_IMAGE[0] NE 0) THEN BEGIN
   PRINT,'CALLING SEQUENCE: MID_UP_IMAGE,IMAGE,DATA,NAXIS,NPIX'
   PRINT,'IMAGE must be a Scalar'
   RETURN
   ENDIF
 
; The regular shutdown procedure for this routine can also handle the shutdown
; required after an error.
 
ON_IOERROR,CLEAN_UP
 
;-------------------------------------------------------------------------------
; Step 1:
; If IMAGE was a string then allocate a LUN and open the file, otherwise assume
; the file is open and use IMAGE as the LUN.
 
IF (INFO_IMAGE[1] EQ 7) THEN BEGIN
   GET_LUN,LUN1
   OPENU,LUN1,IMAGE+'.bdf'
ENDIF ELSE BEGIN
   LUN1 = IMAGE
   ENDELSE
 
;-------------------------------------------------------------------------------
; STEP 2
 
; Read in the directory descriptor information.
 
MID_RD_DIRDSC,LUN1,'NAXIS',NAXIS_IM
IF (NAXIS_IM EQ 0) THEN BEGIN
   PRINT,'ERROR: NAXIS is zero !!!'
   GOTO,CLEAN_UP
   ENDIF
 
MID_RD_DIRDSC,LUN1,'NPIX',NPIX_IM
 
; Truncate NPIX so that its dimension corresponds to NAXIS.
 
IF (NAXIS_IM GT 1) THEN NPIX_IM = NPIX_IM[0:(NAXIS_IM-1)]
 
;-------------------------------------------------------------------------------
; STEP 3
 
; The pixel data is in MIDAS BDF filesection MAINSEG4. The location of MAINSEG4
; (in blocks) is keep in the Frame Control Block (FCB), which is the first block
; of the file.  The MAINSEG list begins at byte 172, consists of I*4 values and
; is 5 values long.  MAINSEG4, thus begins at byte 184 (zero index).
 
F1       = ASSOC(LUN1,BYTARR(512))
MAINSEG4 = LONG(F1[184:187,0],0) - 1    ; subtract one for zero index
; Unix v.2 IDL requires offsets in bytes:
MAINSEG4 = MAINSEG4*512 
;-------------------------------------------------------------------------------
; STEP 4:
 
; Load the MIDAS pixel image starting from the block indicated by MAINSEG4
; using the format indicated by the NAXIS and NPIX descriptors.
 
IF (NAXIS EQ 2) THEN BEGIN
   NP1 = NPIX_IM[0]
   NP2 = NPIX_IM[1]
   F1   = ASSOC(LUN1,FLTARR(NP1,NP2),MAINSEG4)
   ENDIF ELSE BEGIN
   F1   = ASSOC(LUN1,FLTARR(NPIX_IM),MAINSEG4)
   ENDELSE   
DATT = F1[0]
 
;-------------------------------------------------------------------------------
; STEP 5:

; Check for consistency of dimensions before writing data into image.

IF NAXIS NE NAXIS_IM THEN BEGIN
   PRINT, 'ERROR: Dimensions of data array not consistent with image.'
   RETURN
ENDIF

IF (NAXIS GT 2) OR (NAXIS_IM GT 2) THEN BEGIN
   PRINT, 'ERROR: Only 1-d and 2-d cases supported.'
   RETURN
ENDIF

IF NAXIS EQ 1 THEN BEGIN
   IF NPIX LT NPIX_IM THEN BEGIN
      PRINT, 'ERROR: Dimension of data array too large.'
      RETURN
   ENDIF
   IF NPIX GT NPIX_IM THEN BEGIN
      PRINT, 'WARNING: Writing into part (only) of MIDAS image.'
   ENDIF
   DATT[0] = DATA
   F1[0] = DATT
ENDIF

IF NAXIS EQ 2 THEN BEGIN
   IF (NPIX[0] LT NPIX_IM[0]) OR (NPIX[1] LT NPIX_IM[1]) THEN BEGIN
      PRINT, 'ERROR: Dimensions of data array too large.'
      RETURN
   ENDIF
   IF (NPIX[0] GT NPIX_IM[0]) OR (NPIX[1] GT NPIX_IM[1]) THEN BEGIN
      PRINT, 'WARNING: writing into part (only) of MIDAS image.'
   ENDIF
   DATT[0,0] = DATA
   F1[0] = DATT
ENDIF
 
;-------------------------------------------------------------------------------
; STEP 6: clean up.
 
CLEAN_UP:
 
; If the IMAGE input was a string we need to clean up the file action.
 
IF (INFO_IMAGE[1] EQ 7) THEN BEGIN
   FREE_LUN,LUN1
   ENDIF
 
;-------------------------------------------------------------------------------

RETURN
END

