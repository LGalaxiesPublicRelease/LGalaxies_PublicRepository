;-----------------------------------------------------------------------
;+
; NAME:
;   djs_plotlimitbox
;
; PURPOSE:
;   Plot a box that bounds the given limits in X and Y.
;
; CALLING SEQUENCE:
;   djs_plotlimitbox, xrange, yrange
;
; INPUT:
;   xrange:    Range in X
;   yrange:    Range in Y
;
; OUTPUTS:
;
; PROCEDURES CALLED:
;   djs_oplot
;
; REVISION HISTORY:
;   Written by D. Schlegel, 11 Dec 1998, Princeton
;-
;-----------------------------------------------------------------------
pro djs_plotlimitbox, xrange, yrange, _EXTRA=KeywordsForPlot

   ; Set variables to +1 if axis increasing, -1 if decreasing
   xs = 2*(!x.crange[1] GE !x.crange[0]) - 1
   ys = 2*(!y.crange[1] GE !y.crange[0]) - 1

   if (keyword_set(xrange)) then begin
      qLft = xrange[0]*xs GT !x.crange[0]*xs AND xrange[0]*xs LT !x.crange[1]*xs
      qRgt = xrange[1]*xs GT !x.crange[0]*xs AND xrange[1]*xs LT !x.crange[1]*xs
      if (keyword_set(yrange)) then yplot = yrange $
       else yplot = !y.crange

      if (qLft) then begin
         djs_oplot, [xrange[0],xrange[0]], yplot, _EXTRA=KeywordsForPlot
      endif

      if (qRgt) then begin
         djs_oplot, [xrange[1],xrange[1]], yplot, _EXTRA=KeywordsForPlot
      endif
   endif

   if (keyword_set(yrange)) then begin
      qBot = yrange[0]*ys GT !y.crange[0]*ys AND yrange[0]*ys LT !y.crange[1]*ys
      qTop = yrange[1]*ys GT !y.crange[0]*ys AND yrange[1]*ys LT !y.crange[1]*ys
      if (keyword_set(xrange)) then xplot = xrange $
       else xplot = !x.crange

      if (qBot) then begin
         djs_oplot, xplot, [yrange[0],yrange[0]], _EXTRA=KeywordsForPlot
      endif

      if (qTop) then begin
         djs_oplot, xplot, [yrange[1],yrange[1]], _EXTRA=KeywordsForPlot
      endif
   endif

   return
end 
;-----------------------------------------------------------------------
