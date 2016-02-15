;+
; NAME: 
;	DELVARX
; PURPOSE: 
; 	Delete variables for memory management (can call from routines) 
; EXPLANATION:
;	Like intrinsic DELVAR function, but can be used from any calling level
;
; CALLING SEQUENCE:
; 	DELVARX,  a [,b,c,d,e,f,g,h,i,j, /FREE_MEM]
;
; INPUTS: 
;	p0, p1...p9 - variables to delete
;
; OPTIONAL KEYWORD:
;       /FREE_MEM - If set, then free memory associated with pointers 
;                   and objects.   
; RESTRICTIONS: 
;	Can't use recursively due to EXECUTE function
;
; METHOD: 
;	Uses EXECUTE and TEMPORARY function   
;
; REVISION HISTORY:
;	Copied from the Solar library, written by slf, 25-Feb-1993
;	Added to Astronomy Library,  September 1995
;	Converted to IDL V5.0   W. Landsman   September 1997
;       Modified, 26-Mar-2003, Zarro (EER/GSFC) 26-Mar-2003
;       - added FREE_MEM to free pointer/objects
;-

PRO delvarx, p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,free_mem = free_mem

   FOR i = 0, N_PARAMS()-1 DO BEGIN ; for each parameter
      param = STRCOMPRESS("p" + STRING(i),/remove)
;  only delete if defined on inpu (avoids error message)
      exestat = execute("defined=n_elements(" + param + ")" ) 

      IF defined GT 0 THEN BEGIN
         if keyword_set(free_mem) then begin
                  exestat = execute("heap_free," + param)
         ENDIF
         exestat = execute(param + "=0")
         exestat = execute("dvar=temporary(" + param + ")" )
      ENDIF
   ENDFOR
   RETURN
END
