PRO  oploterror, x, y, xerr, yerr, NOHAT=hat, HATLENGTH=hln, ERRTHICK=eth, $
      ERRSTYLE=est, THICK = thick, NOCLIP=noclip, ERRCOLOR = ecol, Nsum = nsum,$
      NSKIP=nskip, LOBAR=lobar, HIBAR=hibar, _EXTRA = pkey, ANONYMOUS_ = Dummy_
;+
; NAME:
;      OPLOTERROR
; PURPOSE:
;      Over-plot data points with accompanying X or Y error bars.
; EXPLANATION:
;      For use instead of PLOTERROR when the plotting system has already been
;      defined. 
;
; CALLING SEQUENCE:
;      oploterror, [ x,]  y, [xerr], yerr,   
;            [ /NOHAT, HATLENGTH= , ERRTHICK =, ERRSTYLE=, ERRCOLOR =, 
;              /LOBAR, /HIBAR, NSKIP = , NSUM = , ... OPLOT keywords ]
; INPUTS:
;      X = array of abcissae, any datatype except string
;      Y = array of Y values, any datatype except string
;      XERR = array of error bar values (along X)
;      YERR = array of error bar values (along Y)
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;      /NOHAT     = if specified and non-zero, the error bars are drawn
;                  without hats.
;      HATLENGTH = the length of the hat lines used to cap the error bars.
;                  Defaults to !D.X_VSIZE / 100).
;      ERRTHICK  = the thickness of the error bar lines.  Defaults to the
;                  THICK plotting keyword.
;      ERRSTYLE  = the line style to use when drawing the error bars.  Uses
;                  the same codes as LINESTYLE.
;      ERRCOLOR =  scalar integer (0 - !D.N_TABLE) specifying the color to
;                  use for the error bars
;      NSKIP = Positive Integer specifying the error bars to be plotted.   
;            For example, if NSKIP = 2 then every other error bar is 
;            plotted; if NSKIP=3 then every third error bar is plotted.   
;            Default is to plot every error bar (NSKIP = 1)
;      NSUM =  Number of points to average over before plotting, default = 
;             !P.NSUM  The errors are also averaged, and then divided by 
;             sqrt(NSUM).   This approximation is meaningful only when the 
;             neighboring error bars have similar sizes.
; 
;      /LOBAR = if specified and non-zero, will draw only the -ERR error bars.
;      /HIBAR = if specified and non-zero, will draw only the +ERR error bars.
;                  If neither LOBAR or HIBAR are set _or_ if both are set,
;                  you will get both error bars.  Just specify one if you
;                  only want one set.
;     Any valid keywords to the OPLOT command (e.g. PSYM, YRANGE) are also 
;     accepted by OPLOTERROR via the _EXTRA facility.
;
; NOTES:
;     If only two parameters are input, they are taken as Y and YERR.  If only
;     three parameters are input, they will be taken as X, Y and YERR, 
;     respectively.
;
; EXAMPLE:
;      Suppose one has X and Y vectors with associated errors XERR and YERR
;      and that a plotting system has already been defined:
;
;       (1) Overplot Y vs. X with both X and Y errors and no lines connecting
;           the points
;                  IDL> oploterror, x, y, xerr, yerr, psym=3
;
;       (2) Like (1) but overplot only the Y errors bars and omits "hats"
;                  IDL> oploterror, x, y, yerr, psym=3, /NOHAT
;
;       (3) Like (2) but suppose one has a positive error vector YERR1, and 
;               a negative error vector YERR2 (asymmetric error bars)
;                  IDL> oploterror, x, y, yerr1, psym=3, /NOHAT,/HIBAR
;                  IDL> oploterror, x, y, yerr2, psym=3, /NOHAT,/LOBAR
;
; PROCEDURE:
;      A plot of X versus Y with error bars drawn from Y - YERR to Y + YERR
;      and optionally from X - XERR to X + XERR is written to the output device
;
; WARNING:
;      This an enhanced version of the procedure OPLOTERR in the standard RSI
;      library.    It was renamed to OPLOTERROR in June 1998 in the IDL 
;      Astronomy library.
;
; MODIFICATION HISTORY:
;      Adapted from the most recent version of PLOTERR.  M. R. Greason,
;            Hughes STX, 11 August 1992.
;      Added COLOR keyword option to error bars W. Landsman   November 1993
;      Add ERRCOLOR, use _EXTRA keyword,           W. Landsman, July 1995
;      Remove spurious call to PLOT_KEYWORDS     W. Landsman, August 1995
;      OPLOT more than 32767 error bars          W. Landsman, Feb 1996
;      Added NSKIP keyword                       W. Landsman, Dec 1996
;      Added HIBAR and LOBAR keywords, M. Buie, Lowell Obs., Feb 1998
;      Rename to OPLOTERROR    W. Landsman    June 1998
;      Converted to IDL V5.0   W. Landsman    June 1998
;      Ignore !P.PSYM when drawing error bars   W. Landsman   Jan 1999
;      Handle NSUM keyword correctly           W. Landsman    Aug 1999
;      Check limits for logarithmic axes       W. Landsman    Nov. 1999
;      Work in the presence of  NAN values     W. Landsman    Dec 2000
;      Improve logic when NSUM or !P.NSUM is set  W. Landsman      Jan 2001
;      Remove NSUM keyword from PLOTS call    W. Landsman      March 2001
;      Only draw error bars with in XRANGE (for speed)  W. Landsman Jan 2002
;      Fix Jan 2002 update to work with log plots  W. Landsman Jun 2002
;-
;                  Check the parameters.
;
 On_error, 2
 np = N_params()
 IF (np LT 2) THEN BEGIN
      print, "OPLOTERR must be called with at least two parameters."
      print, "Syntax: oploterr, [x,] y, [xerr], yerr, [..oplot keywords... "
      print,'     /NOHAT, HATLENGTH = , ERRTHICK=, ERRSTLYE=, ERRCOLOR='
      print,'     /LOBAR, /HIBAR, NSKIP= ]'
      RETURN
 ENDIF

; Error bar keywords (except for HATLENGTH; this one will be taken care of 
; later, when it is time to deal with the error bar hats).

 IF (keyword_set(hat)) THEN hat = 0 ELSE hat = 1
 if not keyword_set(THICK) then thick = !P.THICK
 IF (n_elements(eth) EQ 0) THEN eth = thick
 IF (n_elements(est) EQ 0) THEN est = 0
 IF (n_elements(ecol) EQ 0) THEN ecol = !P.COLOR
 if N_elements( NOCLIP ) EQ 0 THEN noclip = 0
 if not keyword_set(NSKIP) then nskip = 1
 if N_elements(nsum) EQ 0 then nsum = !P.NSUM
 if not keyword_set(lobar) and not keyword_set(hibar) then begin
      lobar=1
      hibar=1
 endif else if keyword_set(lobar) and keyword_set(hibar) then begin
      lobar=1
      hibar=1
 endif else if keyword_set(lobar) then begin
      lobar=1
      hibar=0
 endif else begin
      lobar=0
      hibar=1
 endelse
;
; If no X array has been supplied, create one.  Make sure the rest of the 
; procedure can know which parameter is which.
;
 IF np EQ 2 THEN BEGIN                  ; Only Y and YERR passed.
      yerr = y
      yy = x
      xx = indgen(n_elements(yy))
      xerr = make_array(size=size(xx))

 ENDIF ELSE IF np EQ 3 THEN BEGIN       ; X, Y, and YERR passed.
        yerr = xerr
        yy = y
        xx = x

 ENDIF ELSE BEGIN                        ; X, Y, XERR and YERR passed.
      yy = y
      g = where(finite(xerr))
      xerr[g] = abs(xerr[g])
      xx = x
 ENDELSE

 g = where(finite(yerr))
 yerr[g] = abs(yerr[g])

;
;                  Determine the number of points being plotted.  This
;                  is the size of the smallest of the three arrays
;                  passed to the procedure.  Truncate any overlong arrays.
;

 n = N_elements(xx) < N_elements(yy)

 IF np GT 2 then n = n < N_elements(yerr)   
 IF np EQ 4 then n = n < N_elements(xerr)

 xx = xx[0:n-1]
 yy = yy[0:n-1]
 yerr = yerr[0:n-1]
 IF np EQ 4 then xerr = xerr[0:n-1]

; If NSUM is greater than one, then we need to smooth ourselves (using FREBIN)

 if NSum GT 1 then begin
      n1 = float(n) / nsum
      n  = long(n1)
      xx = frebin(xx, n1)
      yy = frebin(yy, n1)
      yerror = frebin(yerr,n1)/sqrt(nsum)
      if NP EQ 4 then xerror = frebin(xerr,n1)/sqrt(nsum)
  endif else begin
      yerror = yerr
      if NP EQ 4 then xerror = xerr
  endelse

 ylo = yy - yerror*lobar
 yhi = yy + yerror*hibar

 if Np EQ 4 then begin
     xlo = xx - xerror*lobar
     xhi = xx + xerror*hibar
 endif
;
;                  Plot the positions.
;
 if n NE 1 then begin
     oplot, xx, yy, NOCLIP=noclip,THICK = thick,_EXTRA = pkey 
 endif else begin 
     plots, xx, yy, NOCLIP=noclip,THICK = thick,_EXTRA = pkey
 endelse
;
; Plot the error bars.   Compute the hat length in device coordinates
; so that it remains fixed even when doing logarithmic plots.
;
 data_low = convert_coord(xx,ylo,/TO_DEVICE)
 data_hi = convert_coord(xx,yhi,/TO_DEVICE)
 if NP EQ 4 then begin
    x_low = convert_coord(xlo,yy,/TO_DEVICE)
    x_hi = convert_coord(xhi,yy,/TO_DEVICE)
 endif
 ycrange = !Y.CRANGE   &  xcrange = !X.CRANGE
    if !Y.type EQ 1 then ylo = ylo > 10^ycrange[0]
    if (!X.type EQ 1) and (np EQ 4) then xlo = xlo > 10^xcrange[0]
 sv_psym = !P.PSYM & !P.PSYM = 0     ;Turn off !P.PSYM for error bars
; Only draw error bars for X values within XCRANGE
    if !X.TYPE EQ 1 then xcrange = 10^xcrange
    g = where((xx GT xcrange[0]) and (xx LE xcrange[1]), Ng)
    if (Ng GT 0) and (Ng NE n) then begin  
          istart = min(g, max = iend)  
    endif else begin
          istart = 0L & iend = n-1
    endelse
    
 FOR i = istart, iend, Nskip DO BEGIN

    plots, [xx[i],xx[i]], [ylo[i],yhi[i]], LINESTYLE=est,THICK=eth,  $
           NOCLIP = noclip, COLOR = ecol

    ; Plot X-error bars 
    ;
    if np EQ 4 then $
       plots, [xlo[i],xhi[i]],[yy[i],yy[i]],LINESTYLE=est, $
              THICK=eth, COLOR = ecol, NOCLIP = noclip

    IF (hat NE 0) THEN BEGIN
       IF (N_elements(hln) EQ 0) THEN hln = !D.X_VSIZE/100. 
       exx1 = data_low[0,i] - hln/2.
       exx2 = exx1 + hln
       if lobar then $
          plots, [exx1,exx2], [data_low[1,i],data_low[1,i]],COLOR=ecol, $
                 LINESTYLE=est,THICK=eth,/DEVICE, noclip = noclip
       if hibar then $
          plots, [exx1,exx2], [data_hi[1,i],data_hi[1,i]], COLOR = ecol,$
                 LINESTYLE=est,THICK=eth,/DEVICE, noclip = noclip
;                                          
       IF np EQ 4 THEN BEGIN
          IF (N_elements(hln) EQ 0) THEN hln = !D.Y_VSIZE/100.
             eyy1 = x_low[1,i] - hln/2.
             eyy2 = eyy1 + hln
             if lobar then $
                plots, [x_low[0,i],x_low[0,i]], [eyy1,eyy2],COLOR = ecol, $
                       LINESTYLE=est,THICK=eth,/DEVICE, NOCLIP = noclip
             if hibar then $
                plots, [x_hi[0,i],x_hi[0,i]], [eyy1,eyy2],COLOR = ecol, $
                       LINESTYLE=est,THICK=eth,/DEVICE, NOCLIP = noclip
          ENDIF
       ENDIF
    NOPLOT:
ENDFOR
 !P.PSYM = sv_psym 
;
RETURN
END
