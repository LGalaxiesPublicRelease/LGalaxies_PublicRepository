;-----------------------------------------------------------------------
;+
; NAME:
;   djs_oplot
;
; PURPOSE:
;   Modified version of OPLOT.
;
; CALLING SEQUENCE:
;   djs_oplot, [x,] y
;
; INPUT:
;   x:
;   y:
;
; OUTPUTS:
;
; COMMENTS:
;   Allows COLOR, PSYM, and SYMSIZE to be vectors.
;   Also allows COLOR to be string descriptions of eight possible colors.
;   If string descriptors are used, then load a basic 8-color color table.
;
; PROCEDURES CALLED:
;   djs_icolor()
;
; REVISION HISTORY:
;   Written by D. Schlegel, 27 September 1997, Durham
;-
;-----------------------------------------------------------------------
pro djs_oplot, x, y, $
 color=color, psym=psym, symsize=symsize, _EXTRA=KeywordsForPlot

   if (NOT keyword_set(color)) then color = !p.color
   if (NOT keyword_set(psym)) then psym = !p.psym
   if (NOT keyword_set(symsize)) then symsize = 1.0

   ncolor = N_elements(color)
   npsym = N_elements(psym)
   nsize = N_elements(symsize)
   icolor = djs_icolor(color)

   ; If X values don't exist, then create them as PLOT or OPLOT would do
   npt = N_elements(x)
   if (n_elements(y) GT 0) then begin
      xtmp = x
      ytmp = y
   endif else begin
      xtmp = lindgen(npt)
      ytmp = x
   endelse

   if (ncolor LE 1 AND npsym LE 1 AND nsize LE 1) then begin
      oplot, xtmp, ytmp, $
       color=icolor[0], psym=psym[0], symsize=symsize[0], $
       _EXTRA=KeywordsForPlot
   endif else begin
      for ipt=0L, npt-1 do begin
         color1 = icolor[ipt MOD ncolor]
         psym1 = psym[ipt MOD npsym]
         symsize1 = symsize[ipt MOD nsize]
         oplot, [xtmp[ipt]], [ytmp[ipt]], $
          color=color1, psym=psym1, symsize=symsize1, $
          _EXTRA=KeywordsForPlot
      endfor
   endelse

   return
end 
;-----------------------------------------------------------------------
