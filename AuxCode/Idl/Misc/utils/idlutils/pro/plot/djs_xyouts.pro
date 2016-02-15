;-----------------------------------------------------------------------
;+
; NAME:
;   djs_xyouts
;
; PURPOSE:
;   Modified version of XYOUTS
;
; CALLING SEQUENCE:
;   djs_xyouts
;
; INPUT:
;
; OUTPUTS:
;
; COMMENTS:
;   Allows COLOR, and CHARSIZE to be vectors.
;   Also allows COLOR to be string descriptions of eight possible colors.
;   If string descriptors are used, then load a basic 8-color color table.
;
; PROCEDURES CALLED:
;   djs_icolor()
;   TeXtoIDL()
;
; REVISION HISTORY:
;   16-Apr-2000 Written by D. Schlegel, Princeton
;-
;-----------------------------------------------------------------------
pro djs_xyouts, x, y, s, color=color, charsize=charsize, $
 _EXTRA=KeywordsForXyouts

   if (n_elements(color) EQ 0) then color = !p.color
   if (NOT keyword_set(charsize)) then charsize = !p.charsize

   ncolor = N_elements(color)
   nsize = N_elements(charsize)
   icolor = djs_icolor(color)
   npt = n_elements(x)

   if (ncolor LE 1 AND nsize LE 1) then begin
      xyouts, x, y, TeXtoIDL(s), color=icolor, charsize=charsize, $
        _EXTRA=KeywordsForXyouts
   endif else begin
      for ipt=0L, npt-1 do begin
         color1 = icolor[ipt MOD ncolor]
         charsize1 = charsize[ipt MOD nsize]
         xyouts, x, y, TeXtoIDL(s), color=color1, _EXTRA=KeywordsForXyouts
      endfor
   endelse

   return
end 
;-----------------------------------------------------------------------
