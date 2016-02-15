;+
;NAME:
;  nw_rebin_image
;PURPOSE:
;  Divide the dimensions of the image by specified value
;CALLING SEQUENCE:
;  colors= nw_rebin_image(colors,rebin)
;INPUTS:
;  colors      - [NX,NY,3] array containing the R, G, and B images
;  rebin       - factor by which to reduce the size of the output
;                image from the input image; ie, if passed a 200x200
;                image and rebin=2, it will return a 100x100 image; if
;                passed rebin=0.5, it will return a 400x400 image.
;OPTIONAL INPUTS:
;  none
;KEYWORDS:
;  none
;OUTPUTS:
;  The resized image 
;BUGS:
;  Non-robust checking for whether or not to use IDL "rebin".
;DEPENDENCIES:
;
;REVISION HISTORY:
;  11/14/03 written - wherry
;-
FUNCTION nw_rebin_image,colors,rebin1
rebinfactor= 1D0/rebin1
dim = size(colors,/dimensions)
NX_new = round(dim[0]*rebinfactor)
NY_new = round(dim[1]*rebinfactor)
rebinned_colors = fltarr(NX_new,NY_new,3)
if (rebinfactor GT 1) AND (round(rebinfactor) EQ rebinfactor) then begin
    FOR k=0,2 DO $
      rebinned_colors[*,*,k] = rebin(colors[*,*,k],NX_new,NY_new,/sample)
endif else begin
    FOR k=0,2 DO $
      rebinned_colors[*,*,k] = congrid(colors[*,*,k],NX_new,NY_new,/interp)
endelse
RETURN,rebinned_colors
END
