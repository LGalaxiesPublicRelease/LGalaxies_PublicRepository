pro one_ray,xcen,ycen,len,angle,terminus,thick=thick,color=color,nodraw=nodraw
;+
; NAME:
;	ONE_RAY
; PURPOSE:
;	Draw a line with a specified starting point, length, and  angle
;
; CALLING SEQUENCE:
;	one_ray, xcen, ycen, len, angle, terminus, [  THICK=, COLOR =, /NODRAW ]
;
; INPUT PARAMETERS:
;	xcen, ycen = starting point in device coordinates, floating point 
;			scalars
;	len        = length in pixels, device coordinates
;	angle      = angle in degrees counterclockwise from +X direction
;
; OUTPUT PARAMETERS:
;	terminus = two-element vector giving ending point of ray in device
;		coordinates
;
; OPTIONAL KEYWORD INPUT PARAMETERS:
;	thick    usual IDL meaning, default = 1.0
;	color    usual IDL meaning, default = !P.COLOR
;	nodraw   if non-zero, the ray is not actually drawn, but the terminus
;		is still calculated
;
; EXAMPLE:
;	Draw a double thickness line of length 32 pixels from (256,256) 
;	45 degrees counterclockwise from the X axis
;
;	IDL> one_ray, 256, 256, 32, 45 ,term, THICK = 2
;
; PROCEDURE:  straightforward matrix arithmetic
;
; MODIFICATION HISTORY:
;    Written by R. S. Hill, Hughes STX Corp., 20-May-1992.
;    Modified to work correctly for COLOR=0  J.Wm.Parker  HITC   1995 May 25
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_error,2

 if N_params() LT 4 then begin
    print,'Syntax -  one_ray, xcen, ycen, len, angle, terminus, ' + $
               '[  THICK= ,COLOR =, /NODRAW ]'
 endif

 if not keyword_set(thick)     then thick     = 1.0
 if (n_elements(color) eq 0)   then color     = !P.COLOR
 sina = sin(angle/!radeg)
 cosa = cos(angle/!radeg)
 rot_mat = [ [ cosa, sina ], [-sina, cosa ] ]
 terminus =  (rot_mat # [len, 0.0]) + [xcen, ycen]

 if not keyword_set(nodraw) then begin
   plots, [xcen, terminus[0]], [ycen, terminus[1]], $
      color=color,/device,thick=thick
 endif

 return
 end
