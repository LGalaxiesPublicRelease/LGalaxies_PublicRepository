pro tvellipse, rmax, rmin, xc, yc, pos_ang, color, DATA = data, THICK = thick, $
	NPOINTS = npoints, COLOR = thecolor, LINESTYLE = linestyle
;+
; NAME:
;      TVELLIPSE
;
; PURPOSE:
;      Draw an ellipse on the current graphics device.
;
; CALLING SEQUENCE:
;      TVELLIPSE, rmax, rmin, xc, yc, [ pos_ang, color, COLOR= ,/DATA, NPOINTS=
;                                        LINESTYLE=, THICK = 
; INPUTS:
;       RMAX,RMIN - Scalars giving the major and minor axis of the ellipse
; OPTIONAL INPUTS:
;       XC,YC - Scalars giving the position on the TV of the ellipse center
;               If not supplied (or if XC, YC are negative and /DATA is not set), 
;               and an interactive graphics device (e.g. not postscript) is set,
;               then the user will be prompted for X,Y
;       POS_ANG - Position angle of the major axis, measured counter-clockwise
;                 from the X axis.  Default is 0.
;       COLOR - Scalar  giving intensity level to draw ellipse.   The color
;               can be specified either with either this parameter or with the 
;               COLOR keyword.   Default is !P.COLOR
;
; OPTIONAL KEYWORD INPUT:
;        COLOR - Intensity value used to draw the circle, overrides parameter
;               value.  Default = !P.COLOR
;        /DATA - if this keyword is set and non-zero, then the ellipse radii and
;               X,Y position center are interpreted as being in DATA 
;               coordinates.   Note that the data coordinates must have been 
;               previously defined (with a PLOT or CONTOUR call).
;        THICK - Thickness of the drawn ellipse, default = !P.THICK
;        LINESTLYLE - Linestyle used to draw ellipse, default = !P.LINESTYLE
;        NPOINTS - Number of points to connect to draw ellipse, default = 120
;                  Increase this value to improve smoothness
; RESTRICTIONS:
;        TVELLIPSE does not check whether the ellipse is within the boundaries
;        of the window.
;
;        The ellipse is evaluated at NPOINTS (default = 120) points and 
;        connected by straight lines, rather than using the more sophisticated 
;        algorithm used by TVCIRCLE
;
;        TVELLIPSE does not accept normalized coordinates.
;
;        TVELLIPSE is not vectorized; it only draws one ellipse at a time
; EXAMPLE:
;        Draw an ellipse of major axis 50 pixels, minor axis 30 pixels, centered
;        on (250,100), with the major axis inclined 25 degrees counter-clockwise
;        from the X axis.   Use a double thickness line and device coordinates 
;        (default)
;
;	IDL> tvellipse,50,30,250,100,25,thick=2
; NOTES:
;        Note that the position angle for TVELLIPSE (counter-clockwise from the
;        X axis) differs from the astronomical position angle (counter-clockwise
;        from the Y axis). 
;
; REVISION HISTORY:
;        Written  W. Landsman STX          July, 1989            
;        Converted to use with a workstation.  M. Greason, STX, June 1990
;        LINESTYLE keyword, evaluate at 120 points,  W. Landsman HSTX Nov 1995
;        Added NPOINTS keyword, fixed /DATA keyword W. Landsman HSTX Jan 1996
;        Check for reversed /DATA coordinates  P. Mangiafico, W.Landsman May 1996
;        Converted to IDL V5.0   W. Landsman   September 1997
;        Work correctly when X & Y data scales are unequal  December 1998
;        Removed cursor input when -ve coords are entered with /data 
;        keyword set  P. Maxted, Keele, 2002
;-
 On_error,2                              ;Return to caller

 if N_params() lt 2 then begin
   print,'Syntax - TVELLIPSE, rmax, rmin, xc, yc,[ pos_ang, color, COLOR = '
   print,'                         NPOINTS =, LINESTYLE = ,THICK=, /DATA ]'
   return
 endif

 if N_params() lt 4 then $
       cursor, xc, yc, /DEVICE, /NOWAIT      ;Get unroamed,unzoomed coordinates

 if ( (xc LT 0) or (yc LT 0)) and not(keyword_set(data)) then begin
       message,'Position cursor in window ' + strtrim(!D.WINDOW,2) + $
              ' -- then hit mouse button',/INF
       cursor, xc, yc, /DEVICE, /WAIT
         message,'Ellipse is centered at (' + strtrim(xc,2) + ',' + $
		strtrim(yc,2) + ')',/INF
 endif

 if N_params() LT 5 then pos_ang = 0.    ;Default position angle
  if N_Elements(TheColor) EQ 0 then begin
      IF N_Elements( Color ) eq 0 THEN Color = !P.COLOR
  endif else color = TheColor

 if not keyword_set(THICK) then thick = !P.THICK
 if not keyword_set(LINESTYLE) then linestyle = !P.LINESTYLE
 if not keyword_set(NPOINTS) then npoints = 120   ;Number of points to connect
 phi = 2*!pi*(findgen(npoints)/(npoints-1))       ;Divide circle into Npoints
 ang = pos_ang/!RADEG               	          ;Position angle in radians
 cosang = cos(ang)
 sinang = sin(ang)

 x =  rmax*cos(phi)              ;Parameterized equation of ellipse
 y =  rmin*sin(phi)

 xprime = xc + x*cosang - y*sinang   	;Rotate to desired position angle
 yprime = yc + x*sinang + y*cosang

 if keyword_set(data) then $
 plots, xprime, yprime, /DATA, COLOR=color, THICK=thick, LINESTYLE=linestyle $
 else $
 plots, round(xprime), round(yprime), $
        /DEVICE, COLOR=color, THICK=thick, LINESTYLE=linestyle

 return
 end
