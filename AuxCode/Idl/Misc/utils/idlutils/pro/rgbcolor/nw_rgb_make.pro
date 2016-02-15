;+
;NAME:
;  nw_rgb_make
;PURPOSE:
;  Creates JPEG (or TIFF) from images
;CALLING SEQUENCE:
;  nw_rgb_make, Rim, Gim, Bim, [name=, scales=, nonlinearity=, $
;      origin=, rebinfactor=, /saturatetowhite]
;INPUTS:
;  Rim,Gim,Bim - R, G, and B fits file names, or data arrays
;OPTIONAL INPUTS:
;  name        - name of the output jpeg file
;  scales      - (3x1) array to scale the R/G/B
;              - defaults are [1.,1.,1.]
;  nonlinearity- 'b'
;              - b=0 for linear fit, b=Inf for logarithmic
;              - default is 3
;  origin      - (3x1) array containing R0/G0/B0
;              - default is [0,0,0]
;  rebinfactor - integer by which to rebin pixels in the x and y
;                directions; eg, a rebinfactor of 2 halves the number
;                of pixels in each direction and quarters the total
;                number of pixels in the image.
;  quality     - quality input for WRITE_JPEG
;  overlay     - [nx/rebinfactor,ny/rebinfactor,3] image to overlay on
;                the input images
;OPTIONAL KEYWORDS:
;  saturatetowhite  - saturate high-value pixels to white rather than to color
;  tiff        - make tiff instead of jpeg
;  dpitiff     - set TIFF "dots per inch" resolution (only if /tiff set)
;  invert      - ???
;OPTIONAL OUTPUTS:
;  
;EXAMPLE:
;  
;KEYWORDS:
;  none
;OUTPUTS:
;  JPEG (or TIFF)
;DEPENDENCIES:
;  
;BUGS:
;  If the code congridded before making the initial colors matrix, it
;  would use less memory and be faster.
;  
;REVISION HISTORY:
; 12/03/03 written - wherry
;-
PRO nw_rgb_make,Rim,Gim,Bim,name=name,scales=scales,nonlinearity= $
                nonlinearity,origin=origin,rebinfactor=rebinfactor, $
                saturatetowhite=saturatetowhite,quality=quality, $
                overlay=overlay,colors=colors,tiff=tiff,invert=invert, $
                underlay=underlay, dpitiff=dpitiff

;set defaults
IF (keyword_set(tiff)) THEN suffix='tif' ELSE suffix='jpg'
IF (NOT keyword_set(name)) THEN name = 'nw_rgb_make.'+suffix
IF (NOT keyword_set(quality)) THEN quality = 100

;assume Rim,Gim,Bim same type, same size
IF size(rim[0],/tname) eq 'STRING' THEN BEGIN
    R = mrdfits(Rim[0])
    dim = size(R,/dimensions)
    NX = LONG(dim[0])
    NY = LONG(dim[1])
    colors = fltarr(NX,NY,3)
    colors[*,*,0] = temporary(R)
    colors[*,*,1] = mrdfits(Gim[0])
    colors[*,*,2] = mrdfits(Bim[0])
ENDIF ELSE BEGIN
    dim = size(Rim,/dimensions)
    NX = LONG(dim[0])
    NY = LONG(dim[1])
    colors = fltarr(NX,NY,3)
    colors[*,*,0] = Rim
    colors[*,*,1] = Gim
    colors[*,*,2] = Bim
ENDELSE
IF (n_elements(rebinfactor) GT 0) THEN $
  IF (rebinfactor NE 1) THEN $
  colors = nw_rebin_image(temporary(colors),rebinfactor)

colors = nw_scale_rgb(temporary(colors),scales=scales)
splog, 'nw_arcsinh'
colors = nw_arcsinh(temporary(colors),nonlinearity=nonlinearity, /inplace)
IF (NOT keyword_set(saturatetowhite)) THEN BEGIN
    splog, 'nw_cut_to_box'
    colors = nw_cut_to_box(temporary(colors),origin=origin)
ENDIF
IF keyword_set(overlay) THEN colors= (colors > overlay) < 1.0
IF keyword_set(underlay) THEN colors= (colors < (1.-underlay)) > 0.0
splog, 'nw_float_to_byte'
colors = nw_float_to_byte(temporary(colors))
if(keyword_set(invert)) then colors=255-colors

IF keyword_set(tiff) THEN BEGIN
    colors = reverse(temporary(colors),2)
    splog, 'writing tiff'
    WRITE_TIFF,name,planarconfig=2,red=colors[*,*,0],$
      green=colors[*,*,1], blue=colors[*,*,2], xresol=dpitiff, yresol=dpitiff
ENDIF ELSE BEGIN
    splog, 'writing jpeg'
    if(NOT arg_present(colors)) then $
      WRITE_JPEG,name,temporary(colors),TRUE=3,QUALITY=quality $ 
    else $
      WRITE_JPEG,name,(colors),TRUE=3,QUALITY=quality 
ENDELSE
END
