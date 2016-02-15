;-----------------------------------------------------------------------
;+
; NAME:
;   djs_arrow
;
; PURPOSE:
;   Modified version of ARROW to allow a string for the color(s)
;
; CALLING SEQUENCE:
;   djs_arrow
;
; INPUT:
;
; OUTPUTS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Written by D. Schlegel, 11 Dec 1998, Princeton
;-
;-----------------------------------------------------------------------
pro djs_arrow, x0, y0, x1, y1, color=color, _EXTRA=KeywordsForArrow

   icolor = djs_icolor(color)
   arrow, x0, y0, x1, y1, color=icolor, _EXTRA=KeywordsForArrow

   return
end 
;-----------------------------------------------------------------------
