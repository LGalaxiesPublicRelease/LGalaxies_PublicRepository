function abscal,value,header,DEBUG = debug
;+
; NAME:
;       ABSCAL
; PURPOSE:
;       Apply the FITS BZERO and BSCALE keyword values to a data array
;
; CALLING SEQUENCE:
;       RESULT = ABSCAL( Value, Header, /DEBUG)
;
; INPUTS:
;       VALUE -  Any scalar, vector, or array (usually an integer type giving a
;               relative intensity).
;       HEADER - A FITS  header array containing the absolute calibration
;               keyword BSCALE, and optionally BZERO and BUNIT.
;
; OUTPUT:
;       RESULT = BSCALE*VALUE + BZERO, where the BSCALE and BZERO scalars
;               are taken from the FITS header.  
;               If the absolute calibration keywords do not exist, then
;               RESULT = VALUE, and !ERR = -1.
;
; OPTIONAL INPUT KEYWORD:
;       /DEBUG - If DEBUG is set, then ABSCAL will print the
;               calibration units given by the BUNIT keyword.
;
; REVISION HISTORY:
;       Written W. Landsman, STX Corporation     January 1987
;       Use DEBUG keyword instead of !DEBUG      September 1995
;       Converted to IDL V5.0   W. Landsman   September 1997
;-
  On_error,2
 
 if ( N_elements(value) EQ 0 ) then message, $
     'Array (first parameter) to be calibrated is not defined'   
 
 zparcheck, 'ABSCAL', header, 2, 7, 1, 'FITS Image header'
 bscale = sxpar( header, 'BSCALE', Count = N_Bscale)     
 if N_Bscale EQ 0 then begin                 ;Did BSCALE keyword exist?
        print,'ABSCAL - No calibration keywords in FITS header'
        return,value
 endif

 bzero  = sxpar( header, 'BZERO' )          ;Get BZERO value, 0 if not found
 if keyword_set(DEBUG) then begin                 ;Print calibration units          
    bunits = sxpar(header,'BUNIT', Count = N_Bunit)       ;Are calibration units supplied?
   if (N_Bunit GT 0) then message,/INF,'Calibration Units: ' + bunits else $
        message,/INF,'Calibration Units not given'
 endif

 return,value*bscale + bzero        
 
 end
