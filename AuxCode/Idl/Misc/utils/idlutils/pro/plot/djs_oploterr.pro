;-----------------------------------------------------------------------
;+
; NAME:
;   djs_oploterr
;
; PURPOSE:
;   Modified version of OPLOTERR and DJS_OPLOT.
;
;   Allows COLOR, PSYM, and SYMSIZE to be vectors.
;   Also allows COLOR to be string descriptions of eight possible colors.
;   If string descriptors are used, then load a basic 8-color color table.
;
; CALLING SEQUENCE:
;   djs_oploterr, [x,] y, xerr=xerr, yerr=yerr, xlog=xlog, ylog=ylog, $
;    cap=cap, xlen=xlen, ylen=ylen, $
;    color=color, psym=psym, symsize=symsize
;
; INPUT:
;   x:
;   y:
;
; OPTIONAL INPUTS:
;   xerr:   Error in X; or -1 for upper limit arrow, -2 for lower limit arrow
;   yerr:   Error in Y; or -1 for upper limit arrow, -2 for lower limit arrow
;   xlog:   If set, take the logarithm of X and its error
;   ylog:   If set, take the logarithm of Y and its error
;   cap:    If set, place caps on error bars
;   xlen:   Length of upper/lower limit bars in X; default to 6% of plot range
;   ylen:   Length of upper/lower limit bars in Y; default to 6% of plot range
;
; OUTPUTS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Written by D. Schlegel, 5 February 1998, Durham
;-
;-----------------------------------------------------------------------
pro djs_oploterr, x, y, xerr=xerr, yerr=yerr, xlog=xlog, ylog=ylog, $
 cap=cap, xlen=xlen, ylen=ylen, $
 color=color, psym=psym, symsize=symsize, _EXTRA=KeywordsForPlot

   if (NOT keyword_set(color)) then color = !p.color
   if (NOT keyword_set(symsize)) then symsize = 1.0
   if (NOT keyword_set(psym)) then psym = !p.psym
   blah= !p.psym  ; store away psym value for restoration later
   !p.psym= 0

   if (NOT keyword_set(xlen)) then xlen = 0.06 * (!x.crange[1] - !x.crange[0])
   if (NOT keyword_set(ylen)) then ylen = 0.06 * (!y.crange[1] - !y.crange[0])

   ncolor = N_elements(color)
   npsym = N_elements(psym)
   nsize = N_elements(symsize)

   ; If COLOR is a string or array of strings, then convert color names
   ; to integer values
   if (size(color,/tname) EQ 'STRING') then icolor= djs_icolor(color) $
     else icolor= color

   ; If X values don't exist, then create them as PLOT or OPLOT would do
   npt = N_elements(x)
   if (n_elements(y) GT 0) then begin
      xtmp = x
      ytmp = y
   endif else begin
      xtmp = lindgen(npt)
      ytmp = x
   endelse

   ; Find lower and upper limit bars
   if (keyword_set(xerr)) then begin
      if (keyword_set(xlog)) then begin
         xlo = alog10(xtmp-xerr)
         xhi = alog10(xtmp+xerr)
      endif else begin
         xlo = xtmp-xerr
         xhi = xtmp+xerr
      endelse
   endif
   if (keyword_set(yerr)) then begin
      if (keyword_set(ylog)) then begin
         ylo = alog10(ytmp-yerr)
         yhi = alog10(ytmp+yerr)
      endif else begin
         ylo = ytmp-yerr
         yhi = ytmp+yerr
      endelse
   endif

   ; Take the logarithm of X or Y
   if (keyword_set(xlog)) then xtmp = alog10(xtmp)
   if (keyword_set(ylog)) then ytmp = alog10(ytmp)

   for ipt=0L, npt-1 do begin
      color1 = icolor[ipt MOD ncolor]
      psym1 = psym[ipt MOD npsym]
      symsize1 = symsize[ipt MOD nsize]

      ; Plot the point symbol
      oplot, [xtmp[ipt]], [ytmp[ipt]], $
       color=color1, psym=psym1, symsize=symsize1, $
       _EXTRA=KeywordsForPlot

      ; Plot X error bars
      if (keyword_set(xerr)) then begin
         if (xerr[ipt] GT 0) then begin ; Draw limit bars
            djs_oplot, [xlo[ipt],xhi[ipt]], [ytmp[ipt],ytmp[ipt]], psym=0,$
             color=color1, _EXTRA=KeywordsForPlot
         endif else if (xerr[ipt] LT 0) then begin ; Draw upper/lower limit
            xend = xtmp[ipt] - xlen*(xerr[ipt] EQ -1) + xlen*(xerr[ipt] EQ -2)
            arrow, xtmp[ipt], ytmp[ipt], xend, ytmp[ipt], /data, $
             color=color1, _EXTRA=KeywordsForPlot
         endif
      endif

      ; Plot Y error bars
      if (keyword_set(yerr)) then begin
         if (yerr[ipt] GT 0) then begin ; Draw limit bars
            djs_oplot, [xtmp[ipt],xtmp[ipt]], [ylo[ipt],yhi[ipt]], psym=0,$
             color=color1, _EXTRA=KeywordsForPlot
         endif else if (yerr[ipt] LT 0) then begin ; Draw upper/lower limit
            yend = ytmp[ipt] - ylen*(yerr[ipt] EQ -1) + ylen*(yerr[ipt] EQ -2)
;            djs_oplot, [xtmp[ipt],xtmp[ipt]], [ytmp[ipt],yend], $
;             color=color1, _EXTRA=KeywordsForPlot
            arrow, xtmp[ipt], ytmp[ipt], xtmp[ipt], yend, /data, $
             color=color1, _EXTRA=KeywordsForPlot
         endif
      endif

   endfor

   !p.psym= blah  ; restore psym value

   return
end 
;-----------------------------------------------------------------------
