PRO ASTROLIB
;+
; NAME:
;       ASTROLIB
; PURPOSE:
;       Add the non-standard system variables used by the IDL Astronomy Library
; EXPLANATION: 
;       Also defines the environment variable or VMS 
;       logical ASTRO_DATA pointing to the directory containing data files 
;       associated with the IDL Astronomy library (system dependent).
;
; CALLING SEQUENCE:
;       ASTROLIB
;
; INPUTS:
;       None.
;
; OUTPUTS:
;       None.
;
; METHOD:
;       The non-standard system variables !PRIV, !DEBUG, !TEXTUNIT, and 
;       !TEXTOUT are added using DEFSYSV.
;
; REVISION HISTORY:
;       Written, Wayne Landsman, July 1986.
;       Use DEFSYSV instead of ADDSYSVAR           December 1990
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Test for system variable existence before definition    July 2001
;-
  On_error,2   

  defsysv, '!DEBUG', exist = exist
     if not exist then defsysv, '!DEBUG', 0
  defsysv, '!PRIV', exist = exist   
     if not exist then defsysv, '!PRIV', 0    
  defsysv, '!TEXTUNIT', exist = exist
     if not exist then  defsysv, '!TEXTUNIT', 0
  defsysv, '!TEXTOUT', exist = exist 
     if not exist then defsysv, '!TEXTOUT', 1 

; The following code needs to modified for each particular installation

  CASE !VERSION.OS OF
   'MacOS': ;
    'vms' : setlog,'ASTRO_DATA','$1$DUA5:[IDLUSER.DATA]'
    ELSE  : setenv,'ASTRO_DATA=/export/home/ftp/pub/data/'
  ENDCASE

  message,'Astronomy Library system variables have been added',/INF

  return
  end
 
