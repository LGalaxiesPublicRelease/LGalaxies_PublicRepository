;+
;NAME:
;  nw_cut_to_box
;PURPOSE:
;  Limits the pixel values of the image to a 'box', so that the colors
;  do not saturate to white but to a specific color.
;INPUTS:
;  colors      - (NXxNYx3) array that contains the R/G/B images
;OPTIONAL INPUTS:
;  origin      - (3x1) array containing R0/G0/B0
;              - default is [0,0,0]
;KEYWORDS:
;
;OUTPUTS:
;  The color limited image
;BUGS:
;  
;REVISION HISTORY:
;  11/07/03 written - wherry
;  11/12/03 changed default origin - wherry
;-
FUNCTION nw_cut_to_box,colors,origin=origin

if (NOT keyword_set(origin)) then origin = [0.0,0.0,0.0]
if (n_elements(origin) LT 3) then origin = replicate(origin[0],3)

if((n_elements(colors) MOD 3) ne 0) then $
  message,'nw_cut_to_box requires colors to be [nx,ny,3] array'
if(n_elements(colors) eq 3) then $
  colors=reform(colors,1,1,3)
dim = size(colors,/dimensions)
NX = LONG(dim[0])
NY = LONG(dim[1])

pos_dist = 1.-origin
factor = fltarr(NX,NY)
factor = ((colors[*,*,0]/pos_dist[0]) $
               > (colors[*,*,1]/pos_dist[1]) $
               > (colors[*,*,2]/pos_dist[2]))
factor = factor > 1.0

boxed_colors = colors
FOR b=0,2 DO BEGIN
    ; NICK: check this!
    boxed_colors[*,*,b] = origin[b]+colors[*,*,b]/factor
ENDFOR
RETURN,boxed_colors
END
