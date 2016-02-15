PRO spline_smooth,X,Y,Yerr,distance,coefficients,smoothness,xplot,yplot, $
    TEXTOUT=textout,XTITLE=xtitle,YTITLE=ytitle, INTERP=interp,SILENT=silent, $
    PLOT=plot,ERRBAR=errbar
;+
; NAME:
;       SPLINE_SMOOTH
;
; PURPOSE:
;       Compute a cubic smoothing spline to (weighted) data
; EXPLANATION:
;       Construct cubic smoothing spline (or give regression solution) to given
;       data with minimum "roughness" (measured by the energy in the second 
;       derivatives) while restricting the weighted mean square distance
;       of the approximation from the data.  The results may be written to
;       the screen or a file or both and are optionally returned in the 
;       parameters.  The results may be optionally displayed graphically.
;
; CALLING SEQUENCE:
;       SPLINE_SMOOTH,X,Y,Yerr,distance, [coefficients,smoothness,xplot,yplot 
;               [ XTITLE= ,YTITLE=, INTERP=, TEXTOUT=,/SILENT,/PLOT,/ERRBAR]
;
; INPUT PARAMETERS:
;       X - N_POINT element vector containing the data abcissae
;       Y - N_POINT element vector containing the data ordinates
;       Yerr -     estimated uncertainty in ordinates ( positive scalar)
;       distance - upper bound on the weighted mean square distance
;                  of the approximation from the data (non-negative scalar)
;
; OPTIONAL INPUT PARAMETERS
;      xplot -    vector of spline evaluation abcissae
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       TEXTOUT - Controls print output device, defaults to !TEXTOUT
;       
;               textout=1       TERMINAL using /more option
;               textout=2       TERMINAL without /more option
;               textout=3       <program>.prt
;               textout=4       laser.tmp
;               textout=5       user must open file
;               textout = filename (default extension of .prt)
;
; OPTIONAL OUTPUT PARAMETERS:
;       coefficients - N_POINT x 4 element array containing the sequence of
;                      spline coefficients including the smoothed ordinates.
;       smoothness  - N_POINT element vector containing the energy in second 
;                     derivatives of approximated function.
;       yplot       - vector of evaluated spline ordinates.
;
; OPTIONAL OUTPUT KEYWORD PARAMETERS
;       /SILENT     - suppress all printing.
;       /PLOT       - display smooth curve, data ordinates and error bars
;       /ERRBAR     - display error bars
;       XTITLE      - optional title for X-axis
;       YTITLE      - optional title for Y-axis
;       INTERP      - optionally returned interpolated smooth spline
; NOTES:
;       This procedure constructs a smoothing spline according to the method
;       described in "Fundamentals of Image Processing" by A. Jain  [Prentice-
;       Hall : New Jersey 1989].
;       If the distance parameter is sufficiently large a linear regression
;       is performed, otherwise a cubic smoothing spline is constructed.
;
;       This procedure assumes regular sampling and independent identically 
;       distributed normal errors without missing data.  The data are sorted.
;       
;       SPLINE_SMOOTH uses the non-standard system variables !TEXTOUT and
;       !TEXTUNIT.
;       These can be added to one's session using the procedure ASTROLIB.
;
; COMMON BLOCKS:
;       None.
; EXAMPLE:
;       Obtain coefficients of a univariate smoothing spline fitted to data
;       X,Y assuming normally distributed errors Yerr and write the results to
;       a file.
;
;       IDL> SPLINE_SMOOTH, X, Y, Yerr, distance, coefficients, smoothness,
;            t='spline.dat'
;
;       Fit a smoothing spline to observational data.  Suppress all printing 
;       and save the smoothed ordinates in output variables. Display results.
;
;       IDL> SPLINE_SMOOTH, X, Y, Yerr, distance, coefficients, /SILENT, /PLOT
;       
; PROCEDURES CALLED:
;       Procedures TEXTOPEN, TEXTCLOSE, PLOT, PLOTERROR
;
; RESTRCTIONS:
;       This procedure is damn slow and should probably be rewritten using
;       the Cholesky decomposition.
; AUTHOR:
;       Immanuel Freedman (after A. Jain).      December, 1993
; REVISIONS
;       January 12, 1994    I. Freedman (HSTX)  Adjusted formats
;       March   14, 1994    I. Freedman (HSTX)  Improved convergence
;       March   15, 1994    I. Freedman (HSTX)  User-specified interpolates
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Call PLOTERROR instead of PLOTERR  W. Landsman      February 1999
;-

; Return on error
  ON_ERROR,2
; Constants
  FALSE = 0
  TRUE = 1
  lambda = 0.0    ; initial linear regression
  TOLERANCE = 1.0E-3
  MAXIT = 100      ; number of iterations
; Dispatch table
;
IF N_PARAMS() LT 4 THEN BEGIN
 print,'Syntax - SPLINE_SMOOTH,X,Y,Yerr,distance,[coefficients,smoothness,'
 print,'         xplot,yplot,TEXTOUT=,XTITLE=,YTITLE=,INTERP=,/SILENT,/PLOT,
 print,'         /ERRBAR]'
 RETURN
ENDIF

SZ = size(X) & NX = SZ[1]  ;Number of data points
if SZ[0] NE 1 THEN $
BEGIN
 HELP,X
 MESSAGE,'ERROR - Data matrix is not one-dimensional'
ENDIF
SZ = size(Y) & NY = SZ[1]
if SZ[0] NE 1 THEN $
BEGIN
 HELP,Y
 MESSAGE,'ERROR - Data matrix is not one-dimensional'
ENDIF
if NY NE NX THEN $
BEGIN
 HELP,Y
 MESSAGE,'ERROR - Data vector has incorrect number of elements'
ENDIF
IF distance LT 0.0 THEN MESSAGE,'ERROR - negative distance'
IF yerr     LT 0.0 THEN MESSAGE,'ERROR - negative error value'

  N_POINT = NX ; number of data points
  h = (X[N_POINT-1] - X[0])/(N_POINT -1) ;regular sampling without missing data
  IF distance LE TOLERANCE*yerr THEN distance = TOLERANCE*yerr ; interpolation
; Sort data on X
  index = sort(X) & X = X[index] & Y = Y[index]


straight: ; Linear regression [g(x)=a+bx] is solution for large distances 
; Sample averages

  average_x = TOTAL(X)/N_POINT & average_y = TOTAL(Y)/N_POINT
  average_xx = TOTAL(X*X)/N_POINT & average_xy = TOTAL(X*Y)/N_POINT
; Regression solution
  b = (average_xy - average_x*average_y)/(average_xx -average_x*average_x)
  a = average_y - b*average_x
  g = a+b*x & S = TOTAL((Y-g)*(Y-g))  ; sum of squares
  c = 0 & d = 0
  IF S LE distance THEN BEGIN
  MESSAGE,/INFORM,'Distance sufficient for linear regression'
  linear = TRUE
  smoothness = S
  GOTO,display  
  ENDIF
;
   
; The smoothing spline g(x) = a + b.dx + c.dx*dx + d.dx*dx*dx, 0 <= dx <X - X)
;                                                                        i+1 i
linear = FALSE
MESSAGE,/INFORM,'Cubic spline regression solution'
a=fltarr(N_POINT) & b=a & c=a & d=a


; Tridiagonal Toeplitz matrix

Q =fltarr(N_POINT-2,N_POINT-2) ; Jain pp. 296 ff.
Q[lindgen(N_POINT-2)*(N_POINT-1)]   = 4.0
Q[lindgen(N_POINT-3)*(N_POINT-1)+1] = 1.0
Q[lindgen(N_POINT-3)*(N_POINT-1)+N_POINT-2] = 1.0

; Lower triangular Toeplitz matrix

L = fltarr(N_POINT,N_POINT-2)
L[lindgen(N_POINT-2)*(N_POINT+1)]   = 1.0
L[lindgen(N_POINT-2)*(N_POINT+1)+1] = -2.0
L[lindgen(N_POINT-2)*(N_POINT+1)+2] = 1.0

; Define auxiliary matrices
P = (yerr*yerr) *transpose(L)#L
v = transpose(L)#y
; Iterative Newton solution of "smoothness(lambda) = S"
iteration = 0
REPEAT BEGIN
 
ainverse = P + lambda*Q 
SVDC,ainverse,svdw,svdu,svdv          ; Decompose square matrix 
small = WHERE(svdw LT TOLERANCE*MAX(svdw), count)
IF count NE 0 THEN svdw[small] = 1.0 ; Avoid division by zero below threshold
n = N_ELEMENTS(svdw) & svdwp = fltarr(n,n)
svdwp[(n+1)*lindgen(n)] = 1.0/svdw
A = svdv # svdwp # transpose(svdu)

smoothness = transpose(v)#A#P#A#v & smoothness = smoothness[0]
 derivative = 2.0 * transpose(v)#A#Q#A#P#A#v & derivative = derivative[0]
 if derivative EQ 0 THEN MESSAGE,'ERROR - convergence failed. Stop.'
 lambda = lambda + ((smoothness - distance)/derivative) ; correct sign
 OK = (ABS(smoothness - distance) LE TOLERANCE*distance)
iteration = iteration + 1
done = OK OR iteration GE MAXIT
ENDREP UNTIL done
msg = 'ERROR - Convergence failed after ' + string(MAXIT) + ' iterations.'
IF done AND NOT OK THEN BEGIN
 MESSAGE,/INFORM,strcompress(msg) 
 PRINT,'If you want to try again enter a non-negative smoothness...'
 READ,distance
 IF distance GE 0 THEN GOTO, straight
ENDIF
msg = 'Converged after ' + string(iteration) + ' iterations'
MESSAGE,/INFORM,strcompress(msg)


display:         ; evaluate and display results
IF distance LE TOLERANCE*yerr THEN distance = 0.0 ; interpolating splines
IF linear THEN BEGIN
 coefficients = fltarr(2)
 coefficients[0] = a
 coefficients[1] = b
ENDIF ELSE BEGIN

coefficients = fltarr(N_POINT, 4) 
c = lambda*A#v 
a = Y - (yerr*yerr/lambda)*L#c
c = [0.0,c,0.0]
d = (shift(c,-1) - c)/3.0*h & d[0] = 0.0
b = (shift(a,-1) -a)/h -h*c -h*h*d & b[N_POINT-1] = 0

coefficients[*, 0] = a   ; knots at each data point
coefficients[*, 1] = b
coefficients[*, 2] = c
coefficients[*, 3] = d
ENDELSE

; Open output file
IF NOT KEYWORD_SET(TEXTOUT) THEN TEXTOUT = textout
textopen,'SPLINE_SMOOTH',TEXTOUT = textout
printf,!TEXTUNIT,'SPLINE_SMOOTH: '+ SYSTIME()
; print results 
IF NOT KEYWORD_SET(silent) THEN BEGIN
; print formatted coefficients
printf,!TEXTUNIT,' '
 IF linear THEN BEGIN
  sign = " " &  IF b > 0 THEN sign = ' + '
  msg = 'Linear regression: Y = ' + string(FORMAT = '(G10.4)',a) + ' 
  msg = msg + sign + string(FORMAT='(G10.4)',b) + ' X'

 printf,!TEXTUNIT,strcompress(msg)
 ENDIF ELSE BEGIN
 msg = 'Cubic spline regression at knot points: Y = a +bX + cX^2 +d*X^3 '
 printf,!TEXTUNIT,strcompress(msg)
 msg = string(FORMAT = '(6(6X,A1,6X))','X','Y','a','b','c','d')
 printf,!TEXTUNIT,msg
  FOR i=0, N_POINT-1 DO BEGIN
    printf,!TEXTUNIT,X[i],Y[i],a[i],b[i],c[i],d[i]
  ENDFOR
  conf=[N_POINT - sqrt(2.0*N_POINT), N_POINT + sqrt(2.0*N_POINT)]
  msg = 'Confidence interval for smoothness: [' + string(conf[0]) 
  msg = msg + ', ' + string(conf[1]) + ']'
  printf,!TEXTUNIT,strcompress(msg)
 ENDELSE
 msg = 'Distance = '+string(distance) + ', Smoothness = ' + string(smoothness)
 printf,!TEXTUNIT,strcompress(msg)
ENDIF
; Close output file
textclose, TEXTOUT = textout

IF NOT KEYWORD_SET(plot) THEN BEGIN
 IF KEYWORD_SET(interp) THEN MESSAGE,/INF,'INTERP keyword ignored' 
ENDIF ELSE BEGIN

; plot results (piecewise cubic polynomial => GE 4*N_POINT+1)
 SZ = SIZE(xplot)
 sigma = replicate(yerr,N_POINT)
 IF SZ[0] EQ 0 THEN BEGIN
 MESSAGE,/INFORM,'User did not supply evaluation points'
 xplot = (0.25*h)*findgen(4*N_POINT) + x[0] ; regular sampling
ENDIF ELSE MESSAGE,/INFORM,'User supplied evaluation points'

xplot  = xplot[sort(xplot)] ; sort into increasing order
outofrange = WHERE(xplot LT min(X) OR xplot GT max(X),count)
IF count GT 0 THEN BEGIN
 msg = 'WARNING - user supplied evaluation points out of range: '
 FOR loop = 0,N_elements(outofrange) -1 DO BEGIN
  msg = msg + string(xplot[outofrange[loop]])
 ENDFOR
 MESSAGE,/INFORM,strcompress(msg)
ENDIF
inrange = WHERE((xplot GE min(X)) AND (xplot LE max(X)),count)
IF count GT 0 THEN  xplot = xplot[inrange] ; all values in range
yplot  = fltarr(N_elements(xplot))
Xindex = long((xplot-min(X))/h)  ; truncation

IF linear THEN BEGIN
 yplot = b*xplot + a
ENDIF ELSE BEGIN
 tmp = replicate(1.0,4) 
 dx  = xplot - X[xindex]                      ; fractional part
 a   = reform(tmp#a,4*N_POINT,/OVERWRITE)
 b   = reform(tmp#b,4*N_POINT,/OVERWRITE)
 c   = reform(tmp#c,4*N_POINT,/OVERWRITE)
 d   = reform(tmp#d,4*N_POINT,/OVERWRITE)
 a = a[xindex] & b = b[xindex] & c = c[xindex] & d = d[xindex]
 yplot = ((d*dx + c)*dx + b)*dx + a           ; preserve accuracy
ENDELSE
interp=fltarr(N_elements(xplot),2)
interp[*,0] = xplot
interp[*,1] = yplot
; call to PLOTERR
; Label axes as plot of Y versus X, THICK=2

IF NOT KEYWORD_SET(xtitle) THEN XTITLE='X variate'
IF NOT KEYWORD_SET(ytitle) THEN YTITLE = 'Y variate'
if linear THEN BEGIN
 rmsg = 'Linear regression' 
ENDIF ELSE BEGIN
 rmsg = 'Cubic spline regression'
ENDELSE
msg = rmsg +', distance = '+string(FORMAT = '(G10.4)',distance)+', smoothness = '+string(FORMAT = '(G10.4)',smoothness)
msg = strcompress(msg)
TITLE = msg
symbol = 7                        ; x symbol
IF N_POINT GE 20 THEN symbol = 3  ; point
IF KEYWORD_SET(errbar) THEN BEGIN
 PLOTERROR, X, Y, sigma, PSYM = symbol,xtit=xtitle,ytit=ytitle,xsty=2,ysty=2,$
        xthick=2,ythick=2,thick=2,tit=title 
ENDIF ELSE BEGIN
 PLOT, X, Y, PSYM = symbol,xtit=xitle,ytit=ytitle,xsty=2,ysty=2, $
        xthick=2,ythick=2,thick=2,tit=title
ENDELSE
OPLOT,   xplot,yplot              ; smooth curve
ENDELSE
;
RETURN
END
