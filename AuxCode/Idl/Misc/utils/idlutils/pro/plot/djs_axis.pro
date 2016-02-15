;-----------------------------------------------------------------------
;+
; NAME:
;   djs_axis
;
; PURPOSE:
;   Modified version of AXIS
;
; CALLING SEQUENCE:
;   djs_axis
;
; INPUT:
;
; OUTPUTS:
;
; PROCEDURES CALLED:
;   TeXtoIDL()
;
; REVISION HISTORY:
;   Written by D. Schlegel, 21 Jan 1998, Durham
;-
;-----------------------------------------------------------------------
pro djs_axis, xtitle=xtitle, ytitle=ytitle, $
 _EXTRA=KeywordsForAxis

   if (keyword_set(xtitle)) then xtitle_tex = TeXtoIDL(xtitle)
   if (keyword_set(ytitle)) then ytitle_tex = TeXtoIDL(ytitle)

   axis, xtitle=xtitle_tex, ytitle=ytitle_tex, $
    _EXTRA=KeywordsForAxis

   return
end 
;-----------------------------------------------------------------------
