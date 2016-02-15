;+
; NAME:
;   hogg_image_overlay
; PURPOSE:
;   Bitmap PostScript file to make an image overlay.
; CALLING SEQUENCE:
;   overlay= hogg_image_overlay(psfile,naxis1,naxis2)
; EXAMPLES:
; INPUTS:
;   psfile     - input filename
;   naxis1     - width in pixels to make overlay
;   naxis2     - height in pixels
; OPTIONAL INPUTS:
;   factor     - integer factor to use for antialiasing; default 2;
;                set to 1 for no antialiasing
; OUTPUT:
;   overlay    - overlay with plot material added
; COMMENTS:
; BUGS:
;   - Relies on horrifying UNIX bitmapping code.
;   - Makes insecure intermediate PPM file.
; DEPENDENCIES:
;   pstopnm etc
; REVISION HISTORY:
;   2004-02-28  written - Hogg
;-
function hogg_image_overlay, psfile,naxis1,naxis2,factor=factor
if (NOT keyword_set(factor)) then factor=2L
naxis1= round(naxis1)
naxis2= round(naxis2)
cmd= '\pstopnm -ppm -stdout -xborder 0 -yborder 0'+ $
  ' -ysize '+strtrim(string(factor*naxis1),2)+ $
  ' -xsize '+strtrim(string(factor*naxis2),2)+ $
  ' '+psfile+' > '+psfile+'.ppm'
splog, cmd
spawn, cmd
overlay= read_image(psfile+'.ppm')
cmd= '\rm -vf '+psfile+'.ppm'
splog, cmd
spawn, cmd
overlay= float(overlay)/255.0
tmp= rebin(overlay,3,naxis2,naxis1)
overlay= fltarr(naxis1,naxis2,3)
for bb=0,2 do overlay[*,*,bb]= rotate(reform(tmp[bb,*,*],naxis2,naxis1),6)
help, overlay
return, overlay
end
