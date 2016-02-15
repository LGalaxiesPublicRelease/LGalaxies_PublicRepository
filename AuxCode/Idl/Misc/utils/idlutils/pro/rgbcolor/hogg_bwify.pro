;+
; NAME:
;   hogg_bwify
; PURPOSE:
;   make b/w JPG image from color JPG image
; INPUTS:
;   infile       - file name of color image
;   outfile      - file name for b/w image
; OPTIONAL INPUTS:
;   rebinfactor  - make outfile smaller by factor rebinfactor
; BUGS:
;   - Just uses the r image; it should have various options for how to
;     make the b/w image from the color one.
;   - Always inverts; this should be an option.
; REVISION HISTORY:
;   2004-01-07  commented - Hogg
;-
pro hogg_bwify, infile,outfile,rebinfactor=rebinfactor
quality= 90
read_jpeg, infile,colors,true=3
colors[*,*,1]= 255-colors[*,*,1]
colors[*,*,0]= colors[*,*,1]
colors[*,*,2]= colors[*,*,1]
if keyword_set(rebinfactor) then begin
    colors= nw_rebin_image(colors,rebinfactor)
    colors= byte(colors)
endif
WRITE_JPEG, outfile,temporary(colors),TRUE=3,QUALITY=quality
return
end
