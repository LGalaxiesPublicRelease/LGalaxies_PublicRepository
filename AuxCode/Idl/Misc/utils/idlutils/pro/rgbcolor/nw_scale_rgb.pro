;+
;NAME:
;  nw_scale_rgb
;PURPOSE:
;  mulitiplies the RGB image by their respective scales
;CALLING SEQUENCE:
;  nw_scale_rgb, colors, [scales=]
;INPUTS:
;  colors      - (NXxNYx3) array containing the R, G, and B
;OPTIONAL INPUTS:
;  scales      - (3x1) array to scale the R/G/B
;              - defaults are [1.,1.,1.]
;KEYWORDS:
;  none
;OUTPUTS:
;  The RGB image 
;COMMENTS:
;  The input image must be background-subtracted (ie, have zero background).
;BUGS:
;  
;DEPENDENCIES:
;
;REVISION HISTORY:
;  11/07/03 written - wherry
;-
FUNCTION nw_scale_rgb,colors,scales=scales

;set default scales
IF NOT keyword_set(scales) THEN scales = [1.,1.,1.]
;get dimensions
tmp= size(colors,/dimensions)
NX= LONG(tmp[0])
NY= LONG(tmp[1])

scaled_colors= colors
FOR k=0,2 DO scaled_colors[*,*,k] = colors[*,*,k]*scales[k]

RETURN,scaled_colors
END
