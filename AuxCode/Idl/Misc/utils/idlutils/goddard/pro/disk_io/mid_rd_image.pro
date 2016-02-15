PRO MID_RD_IMAGE,IMAGE,DATA,NAXIS,NPIX
;-------------------------------------------------------------------------------
;+
; NAME:
;	MID_RD_IMAGE
;
; PURPOSE:
;	Get a pixel matrix and some support information from a MIDAS file.
;
; CALLING SEQUENCE:
;	MID_RD_IMAGE,IMAGE,DATA,NAXIS,NPIX
;
;INPUTS:
;	IMAGE = Filename or Logical Unit Number.
;	* If a filename is given, the file will be opened and closed using a
;		local LUN.  The filename is that of the MIDAS image, without
;		extension (.BDF will is assumed) or version number (latest 
;		version is assumed).
;	* If a LUN is given, the file associated with that LUN will be used.
;
; OUTPUTS:
;	NAXIS = Number of dimensions in MIDAS image. I*4 values.
;	NPIX  = Array containing the dimensions. I*4 values.
;	DATA  = Array containing the MIDAS image. Dimensions are defined by
;		NAXIS and NPIX. R*4 values.
;
; ALGORITHM:
;	0) Check inputs and set error handling.
;	1) Open file for access using the access method indicated by the type of
;	the input parameter IMAGE.
;	2) Get the descriptors of the IMAGE.
;		a) NAXIS set the output parameter NAXIS to this value.
;		b) NPIX  set the output parameter NPIX to this value.
;	3) Locate the pixel data start block
;	4) Load pixel data into output parameter DATA
;	5) Terminate file access as is proper for the type of parameter IMAGE.
;
; RESTRICTIONS:
;	1) There must be only one FCB and it must be in VBN 1.
;	2) The LDBs must begin in VBN 2.
;	3) All descriptors must be in the first LDB.
;	4) Currently only works for real data, does not check to see if this is
;		true or not.
;	5) Midas file extensions (.bdf, .tbl) assumed lower case.
;
;AUTHORS:
;	FM   - F. Murtagh, ST-ECF
;	SAV  - Stephen A. Voels, USM/DAN
;
;MODIFICATION HISTORY:
;	OCT 1988 FM   Initial programing and decoding of MIDAS files.
;	FEB 1989 SAV  Name and calling sequence change.
;		General reprograming for effeciency and modularity.
;		Additional parameter checking.
;	MAY 1989 FM   Minor change for case of 1-d images. 
;	FEB 1991 FM   Conversion to V. 2 IDL, Unix.
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-------------------------------------------------------------------------------
; STEP 0:
 
; Set output defaults
 
NAXIS =  0L
NPIX  = [0L]
DATA  = [0.0]
 
; Verify the number of parameters and the types of the inputs.
 
N = N_PARAMS(0)
IF (N NE 4) THEN BEGIN
   PRINT,'CALLING SEQUENCE: MID_RD_IMAGE,IMAGE,DATA,NAXIS,NPIX'
   RETURN
   ENDIF
 
INFO_IMAGE = SIZE(IMAGE)
IF (INFO_IMAGE[0] NE 0) THEN BEGIN
   PRINT,'CALLING SEQUENCE: MID_RD_IMAGE,IMAGE,DATA,NAXIS,NPIX'
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
   OPENR,LUN1,IMAGE+'.bdf'
ENDIF ELSE BEGIN
   LUN1 = IMAGE
   ENDELSE
 
;-------------------------------------------------------------------------------
; STEP 2
 
; Read in the directory descriptor information.
 
MID_RD_DIRDSC,LUN1,'NAXIS',NAXIS
IF (NAXIS EQ 0) THEN BEGIN
   PRINT,'ERROR: NAXIS is zero !!!'
   GOTO,CLEAN_UP
   ENDIF
 
MID_RD_DIRDSC,LUN1,'NPIX',NPIX
 
; Truncate NPIX so that its dimension corresponds to NAXIS.
 
IF (NAXIS GT 1) THEN NPIX = NPIX[0:(NAXIS-1)]
 
;-------------------------------------------------------------------------------
; STEP 3
 
; The pixel data is in MIDAS BDF filesection MAINSEG4. The location of MAINSEG4
; (in blocks) is kept in the Frame Control Block (FCB), which is the first block
; of the file.  The MAINSEG list begins at byte 172, consists of I*4 values and
; is 5 values long.  MAINSEG4, thus begins at byte 184 (zero index).
 
F1       = ASSOC(LUN1,BYTARR(512))
MAINSEG4 = LONG(F1[184:187,0],0) - 1    ; subtract one for zero index
; Unix version 2 of IDL requires offset in bytes:
MAINSEG4 = MAINSEG4*512 
;-------------------------------------------------------------------------------
; STEP 4:
 
; Load the MIDAS pixel image starting from the block indicated by MAINSEG4
; using the format indicated by the NAXIS and NPIX descriptors.
IF (NAXIS EQ 2) THEN BEGIN
   NP1 = NPIX[0]
   NP2 = NPIX[1]
   F1   = ASSOC(LUN1,FLTARR(NP1,NP2),MAINSEG4)
   ENDIF ELSE BEGIN
   F1   = ASSOC(LUN1,FLTARR(NPIX),MAINSEG4)
   ENDELSE   
DATA = F1[0]
 
;-------------------------------------------------------------------------------
; STEP 5: clean up.
 
CLEAN_UP:
 
; If the IMAGE input was a string we need to clean up the file action.
 
IF (INFO_IMAGE[1] EQ 7) THEN BEGIN
   FREE_LUN,LUN1
   ENDIF
 
;-------------------------------------------------------------------------------

RETURN
END
