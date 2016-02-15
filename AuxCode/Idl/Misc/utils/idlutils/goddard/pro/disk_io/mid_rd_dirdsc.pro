PRO MID_RD_DIRDSC,IMAGE,DSCNAME,DSCVALUE
;-------------------------------------------------------------------------------
;+
; NAME:
;	MID_RD_DISDSC
; PURPOSE:
;	Get a MIDAS directory descriptor from a MIDAS BDF. 
; EXPLANATION:
;	Note: PORTABLE MIDAS.
;
; CALLING SEQUENCE:
;	MID_RD_DIRDSC,IMAGE,DSCNAME,DSCVALUE
;
; INPUTS:
;	IMAGE = Filename or Logical Unit Number.
;		* If a filename is given, the file will be opened and closed 
;		using a local LUN.  The filename is that of the MIDAS image, 
;		without extension (.BDF will is assumed) or version number 
;		(latest version is assumed).
;		* If a LUN is given, the file associated with that LUN will be
;		 used.
;	DSCNAME = Name of the Descriptor wanted.
;
; OUTPUTS:
;	DSCVALUE = Value of the directory descriptor wanted.
;
; ALGORITHM:
;	0) Check inputs and set error handling.
;	1) Open file for access using the access method indicated by the type of
;		the input parameter IMAGE.
;	2) Find the descriptor by name (string type).
;	3) Decode the 30 byte descriptor block.
;	4) Use #3 to find descriptor data area.
;	5) Convert descriptor data as indicated by descriptor block information.
;	6) Terminate file access as is proper for the type of parameter IMAGE.
;
; RESTRICTIONS:
;	   1) There must be only one FCB and it must be in VBN 1.
;	2) The LDBs must begin in VBN 2.
;	3) All descriptors must be in the first LDB.
;	4) Note: .bdf and .tbl extensions assumed lower case.
;
; AUTHORS:
;   SAV  - Stephen A. Voels, USM/DAN
;
; MODIFICATION HISTORY:
;	FEB 1989 SAV  Initial programming.
;	MAY 1989 FM   Change of name of routine; some minor alterations.
;	AUG 1990 FM   Updates for Portable MIDAS (noted as comments below).
;	FEB 1991 FM   Conversion to V. 2 of IDL, Unix.
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-------------------------------------------------------------------------------
 
; Set output defaults
 
DSCVALUE = ''
 
; Verify the number of parameters and the types of the inputs.
 
N = N_PARAMS(0)
IF (N NE 3) THEN BEGIN
   PRINT,'CALLING SEQUENCE: MID_RD_DIRDSC,IMAGE,DSCNAME,DSCVALUE'
   RETURN
   ENDIF
 
INFO_IMAGE = SIZE(IMAGE)
IF (INFO_IMAGE[0] NE 0) THEN BEGIN
   PRINT,'CALLING SEQUENCE: MID_RD_DIRDSC,IMAGE,DSCNAME,DSCVALUE'
   PRINT,'IMAGE must be a number or string'
   RETURN
   ENDIF
 
INFO_DSCNAME = SIZE(DSCNAME)
IF (INFO_DSCNAME[0] NE 0) OR (INFO_DSCNAME[1] NE 7) THEN BEGIN
   PRINT,'CALLING SEQUENCE: MID_RD_DIRDSC,IMAGE,DSCNAME,DSCVALUE'
   PRINT,'DSCNAME must be a string'
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
 
; The descriptor information begins in the second block of the file and is
; stored within MIDAS Local Descriptor Blocks (LDBs) using DeSCriptor
; DIRectory (DSCDIR) entries.  A LDB is 2048 bytes long and begins with four
; I*4 integers (16 bytes) of housekeeping information: the 512 byte block
; number (one indexed) (Virtual Block Number, VBN), the number of bytes used
; for descriptors, the VBN of the next LDB containing descriptor information,
; and starting index within the VBN.  Each DSCDIR entry is 30 bytes long and
; the entry list begins after the housekeeping information.
 
; Pad a temporary variable with spaces so that we do not have to trim all the
; strings we want to compare.
; Update by FM (Aug. 90): DSCNAME was padded with 15 blanks, and the last
; parameter to STRMID was 15.  Apparently some junk was contained towards the
; end of these spaces.  Descriptors are most often 5 chars. long, so allowing
; for 7 should be adequate.
 
DSCNAME_TEMP = STRMID(DSCNAME + '       ',0,7)
 
; Setup access to the first LDB.
 
F1   = ASSOC(LUN1,BYTARR(2048),512)
LDB1 = F1[0]
END_DSC = LONG(LDB1,4) + 16
 
;-------------------------------------------------------------------------------
; STEP
 
; Search for the proper entry.  The search starts after the housekeeping data
; (J1 = 16).
 
J1 = 16
; Update by FM (Aug. 90): J1:J1+6 part of string searched.  See comment re.
; length of descriptors above.
 WHILE FIX(STRING(LDB1[J1:J1+6]) NE DSCNAME_TEMP) AND (J1 LT END_DSC) DO $
    J1 = J1 + 30

IF (J1 GE END_DSC) THEN BEGIN
   PRINT,'ERROR: Descriptor '+DSCNAME+' was not found'
   GOTO,CLEAN_UP
   ENDIF
 
; Decode the 30 byte descriptor.
; The offset within the DIRDSCs (one indexed) is converted to a byte offset
; within the LDB (zero indexed) by taking into account the 4 bytes per I*4
; unit of the DIRDSCx and 16 bytes of housekeeping data at the begining of
; the LDB.
 
DSC_NAME   = STRING(LDB1[J1:J1+14])
; Update by FM (Aug90): TYPE was 5 chars in pre-portable Midas. Now Klaus has
; changed this to 1 char.  
DSC_TYPE   = STRING(LDB1[J1+15:J1+15])
DSC_BPE    =    FIX(LDB1[J1+20:J1+21],0)
DSC_NUME   =    FIX(LDB1[J1+22:J1+23],0)
DSC_BLK    =   LONG(LDB1[J1+24:J1+27],0)
DSC_OFFSET =    FIX(LDB1[J1+28:J1+29],0)*4 + 12
 
;-------------------------------------------------------------------------------
 
; If the descriptor value is located within the first LDB then it is not
; necessary to read in any data, otherwise we read in the correct LDB.
 
IF (DSC_BLK EQ LONG(LDB1,0)) THEN BEGIN
   DSC_BYTES = LDB1[DSC_OFFSET:DSC_OFFSET+DSC_BPE*DSC_NUME-1]
ENDIF ELSE BEGIN
   LDBX = F1[(DSC_BLK-2)/4]
   DSC_BYTES = LDBX[DSC_OFFSET:DSC_OFFSET+DSC_BPE*DSC_NUME-1]
   ENDELSE
 
;-------------------------------------------------------------------------------
 
; Convert the descriptor data section into the form indicated by the descriptor
; information.
 
CASE DSC_TYPE OF
; Update by FM (Aug90): strings in the following have been shortened to 1 char
; (see remark above).  Also a 'D' case has been added.  This has not been 
; tested!
   'C': DSCVALUE = STRING(DSC_BYTES)
   'I': IF (DSC_NUME EQ 1) THEN DSCVALUE =   LONG(DSC_BYTES,0) $
                               ELSE DSCVALUE =   LONG(DSC_BYTES,0,DSC_NUME)
   'R': IF (DSC_NUME EQ 1) THEN DSCVALUE =  FLOAT(DSC_BYTES,0) $
                               ELSE DSCVALUE =  FLOAT(DSC_BYTES,0,DSC_NUME)
   'D': IF (DSC_NUME EQ 1) THEN DSCVALUE =  DOUBLE(DSC_BYTES,0) $
                               ELSE DSCVALUE =  DOUBLE(DSC_BYTES,0,DSC_NUME)
   ELSE: BEGIN
      PRINT,'ERROR: Unknown type of Descriptor Data'
      GOTO,CLEAN_UP
      END
   ENDCASE
 
 
;-------------------------------------------------------------------------------
 
CLEAN_UP:
 
; If the IMAGE input was a string we need to clean up the file action.
 
IF (INFO_IMAGE[1] EQ 7) THEN BEGIN
   FREE_LUN,LUN1
   ENDIF
 
;-------------------------------------------------------------------------------
 
RETURN
END
 
