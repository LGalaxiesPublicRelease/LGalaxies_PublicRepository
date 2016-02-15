FUNCTION FXMOVE, UNIT, N_EXT, SILENT = Silent

;+
; NAME:
;     FXMOVE
; PURPOSE:
;     Skip a specified number of extensions in a FITS file
;
; CALLING SEQUENCE:
;     STATUS=FXMOVE(UNIT, N_EXT, /Silent)
;
; INPUT PARAMETERS:
;     UNIT     = An open unit descriptor for a FITS data stream.
;     N_EXT   = Number of extensions to skip.
;
; OPTIONAL INPUT PARAMETER:
;     /SILENT - If set, then any messages about invalid characters in the 
;               FITS file are suppressed.
; RETURNS:
;     0 if successful.
;    -1 if an error is encountered.
;
; COMMON BLOCKS:
;      None.
; SIDE EFFECTS:
;      Repositions the file pointer.
; PROCEDURE:
;      Each FITS header is read in and parsed, and the file pointer is moved
;      to where the next FITS extension header until the desired
;      extension is reached.
; PROCEDURE CALLS:
;      FXPAR(), MRD_HREAD, MRD_SKIP
; MODIFICATION HISTORY:
;      Extracted from FXPOSIT 8-March-2000 by T. McGlynn
;      Added /SILENT keyword  14-Dec-2000 by W. Landsman
;      Save time by not reading the full header  W. Landsman   Feb. 2003
;-
  
        FOR EXT = 1, N_EXT DO BEGIN
               
;
;  Read the next header, and get the number of bytes taken up by the data.
; 
                IF EOF(UNIT) THEN RETURN, -1
           
                ; Can't use FXHREAD to read from pipe, since it uses
                ; POINT_LUN.  So we read this in ourselves using mrd_hread

                MRD_HREAD, UNIT, HEADER, STATUS, SILENT = Silent, /FIRSTBLOCK
                IF STATUS LT 0 THEN RETURN, -1
                        
                ; Get parameters that determine size of data
                ; region.
                
                BITPIX = FXPAR(HEADER,'BITPIX')
                NAXIS  = FXPAR(HEADER,'NAXIS')
                GCOUNT = FXPAR(HEADER,'GCOUNT') 
                IF GCOUNT EQ 0 THEN GCOUNT = 1
                PCOUNT = FXPAR(HEADER,'PCOUNT')
                
                IF NAXIS GT 0 THEN BEGIN 
                        DIMS = FXPAR(HEADER,'NAXIS*')           ;Read dimensions
                        NDATA = DIMS[0]
                        IF NAXIS GT 1 THEN FOR I=2,NAXIS DO NDATA = NDATA*DIMS[I-1]
                        
                ENDIF ELSE NDATA = 0
                
                NBYTES = (ABS(BITPIX) / 8) * GCOUNT * (PCOUNT + NDATA)
;
;  Move to the next extension header in the file.
;
                NREC = LONG((NBYTES + 2879) / 2880)
                
                MRD_SKIP, UNIT, NREC*2880L 

        ENDFOR
        
        RETURN, 0
        
END
