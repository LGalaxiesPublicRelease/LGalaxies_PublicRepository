;-----------------------------------------------------------------------
;+
; NAME:
;   djs_contourpts
;
; PURPOSE:
;   Make a contour plot from point data, drawing contours only where the
;   point density is high.
;
; CALLING SEQUENCE:
;   djs_contourpts
;
; INPUT:
;   xpt:
;   ypt:
;
; OPTIONAL KEYWORDS:
;   bin1:
;   bin2:
;   overplot:  If set, then use current plot limits and overplot.
;   nlevels:
;   levels:
;   loglevels: If set, then select NLEVEL (or 6) log-spaced levels
;   nopoints:  If set, then do not plot any point data (only contours).
;   psym:      Keyword for plotting point data; default to 3
;   color:     Keyword for plotting point data
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   level0:    Lowest contour level
;   ilow:      Indices for points outside the lowest contour level.
;
; PROCEDURES CALLED:
;   djs_icolor()
;   djs_oplot
;
; REVISION HISTORY:
;   Written by D. Schlegel, 9 Dec 1998, Princeton
;-
;-----------------------------------------------------------------------
pro djs_contourpts, xpt, ypt, bin1=bin1, bin2=bin2, $
 nlevels=nlevels, levels=levels, loglevels=loglevels, $
 xrange=xrange, yrange=yrange, overplot=overplot, $
 psym=psym, color=color, nopoints=nopoints, ilow=ilow, level0=level0, $
 _EXTRA=KeywordsForPlot

   ; Determine plot ranges
   if (keyword_set(overplot)) then begin
      xrange = !x.crange
      yrange = !y.crange
   endif
   if (NOT keyword_set(xrange)) then xrange = [ min(xpt), max(xpt) ]
   if (NOT keyword_set(yrange)) then yrange = [ min(ypt), max(ypt) ]

   if (NOT keyword_set(bin1)) then bin1 = 1
   if (NOT keyword_set(bin2)) then bin2 = 1

   ; Create the image which is a density map of the point data
   if (xrange[1] GT xrange[0]) then xbin = bin1 $
    else xbin = -bin1
   ix = (xpt-xrange[0]) / xbin
   ixmax = fix( (xrange[1]-xrange[0])/xbin )
   xvec = (indgen(ixmax+1)+0.5) * xbin + xrange[0]
   if (yrange[1] GT yrange[0]) then ybin = bin2 $
    else ybin = -bin2
   iy = (ypt-yrange[0]) / ybin
   iymax = fix( (yrange[1]-yrange[0])/ybin )
   yvec = (indgen(iymax+1)+0.5) * ybin + yrange[0]
   ivalid = where( (ix GE 0) AND (ix LT ixmax) AND (iy GE 0) AND (iy LT iymax) )
   image = hist_2d( ix[ivalid], iy[ivalid], max1=ixmax, max2=iymax )

   ; Determine the NLEVEL, LEVELS if not set
   if (NOT keyword_set(nlevels) AND NOT keyword_set(levels)) then nlevels = 6
   if (NOT keyword_set(levels)) then begin
      if (keyword_set(loglevels)) then begin
         levels = 10^( (findgen(nlevels)+1) * alog10(max(image)) / nlevels )
      endif else begin
         levels = (findgen(nlevels)+1) * max(image) / nlevels
      endelse
   endif

   ; Draw the contour plot
   if (NOT keyword_set(overplot)) then begin
      plot, [0], [0], /nodata, $
       xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, $
       _EXTRA=KeywordsForPlot
   endif
   icolor = djs_icolor(color)
   contour, image, xvec, yvec, /overplot, $
    nlevels=nlevels, levels=levels, color=icolor, c_colors=icolor

   ; Determine the indices for points outside the lowest contour
   level0 = levels[0]
   ilow = where( image(ix[ivalid], iy[ivalid]) LT level0, Nlow )
   if (Nlow GT 0) then begin
      ilow = ivalid[ilow]

      ; Overplot the points outside of the lowest contour
      if (NOT keyword_set(nopoints)) then begin
         if (NOT keyword_set(psym)) then psym=3
;         plot, [0], [0], /nodata, /noerase, $
;          xrange=xrange, yrange=yrange, xstyle=5, ystyle=5, $
;          _EXTRA=KeywordsForPlot
         djs_oplot, xpt[ilow], ypt[ilow], psym=psym, color=color, $
          _EXTRA=KeywordsForPlot
      endif
   endif

   return
end 
;-----------------------------------------------------------------------
