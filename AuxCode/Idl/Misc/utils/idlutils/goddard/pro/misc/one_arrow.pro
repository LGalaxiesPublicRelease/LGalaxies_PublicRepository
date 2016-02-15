pro one_arrow,xcen,ycen,angle,label, $
              charsize=charsize,thick=thick,color=color, $
              arrowsize=arrowsize,font = font
;+
; NAME:
;	ONE_ARROW
; PURPOSE:
;	Draws an arrow labeled with a single character on the current device
; EXPLANATION:
;	ONE_ARROW is called, for example, by ARROWS to create a
;	"weathervane" showing the N-E orientation of an image.
;
; CALLING SEQUENCE:
;	one_arrow, xcen, ycen, angle, label, CHARSIZE = , THICK = , COLOR = 
;			ARROWSIZE=, FONT =  ]
; INPUT PARAMETERS:
;    xcen, ycen = starting point of arrow in device coordinates, floating
;			point scalars,
;    angle      = angle of arrow in degrees counterclockwise from +X direction
;    label      = single-character label (may be blank)
;
; OUTPUT PARAMETERS:  none
;
; OPTIONAL INPUT PARAMETERS:
;      	CHARSIZE   = usual IDL meaning, default = 2.0
;	THICK      = usual IDL meaning, default = 2.0
;	COLOR      = usual IDL meaning, default = !P.COLOR
;	ARROWSIZE  = 3-element vector defining appearance of arrow.
;		Default = [30.0, 9.0, 35.0], meaning arrow is 30 pixels
;		long; arrowhead lines 9 pixels long and inclined 35
;		degrees from arrow shaft.
;		If you try to use a non-TV device, you will probably
;		want to change this.
;	FONT - IDL vector font number to use (1-20).   For example, to write
;		the 'N' and 'E' characters in complex script, set font=13
; EXAMPLE:
;	Draw an triple size arrow emanating from the point (212,224)
;	and labeled with the character 'S'
;
;	IDL> one_arrow,212,224,270,'S',charsize=3
; PROCEDURE:  
;	Calls one_ray to vector-draw arrow.
; MODIFICATION HISTORY:
;	Written by R. S. Hill, Hughes STX Corp., 20-May-1992.
;	Added font keyword, W.B. Landsman Hughes STX Corp. April 1995
;	Modified to work correctly for COLOR=0  J.Wm.Parker, HITC   1995 May 25
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_error,2

 if N_params() LT 4 then begin
      print,'Syntax - one_arrow, xcen, ycen, angle, label, CHARSIZE = , ' 
      print,'                         THICK = , COLOR =, ARROWSIZE = ]'
      return
 endif 

 if not keyword_set(charsize)  then charsize  = 2.0
 if not keyword_set(thick)     then thick     = 2.0
 if (n_elements(color) eq 0)   then color     = !P.COLOR
 if (n_elements(arrowsize) ge 1) and (n_elements(arrowsize) ne 3) then begin
   print,'Error in ONE_ARROW:  returning to main level.'
   print,'Arrowsize is [length, head_length, head_angle]'
   print,'Defaults are [30.0,9.0,35.0]'
   retall
 endif

 if n_elements(arrowsize) eq 0 then arrowsize=[30.0,9.0,35.0]
 label = strmid(strtrim(label,2),0,1)
 if keyword_set(font) then label = '!' + strtrim(font,2) + label + '!X '
 len       = arrowsize[0]
 headlen   = arrowsize[1]
 headangle = arrowsize[2]
 baseline  = (!d.y_ch_size+!d.x_ch_size)/2.0
 char_cen_offset  = baseline*charsize
 char_orig_len    = char_cen_offset/2.0
 char_orig_angle  = 225.0

;  Draw shaft of arrow
one_ray,xcen,ycen,len,angle,terminus,thick=thick,color=color

;  Draw head of arrow
one_ray,terminus[0],terminus[1],headlen,(angle+180.0+headangle),t2, $
   thick=thick,color=color
one_ray,terminus[0],terminus[1],headlen,(angle+180.0-headangle),t2, $
   thick=thick,color=color

;  Draw label
one_ray,xcen,ycen,len+char_cen_offset,angle,terminus,/nodraw
one_ray,terminus[0],terminus[1],char_orig_len,char_orig_angle,char_orig,/nodraw
xyouts, char_orig[0], char_orig[1], label, /device, $
   charthick=thick, color=color, charsize=charsize

 return
 end
