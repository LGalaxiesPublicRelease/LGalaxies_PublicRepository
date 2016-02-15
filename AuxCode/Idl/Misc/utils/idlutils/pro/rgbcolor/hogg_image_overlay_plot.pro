;+
; NAME:
;   hogg_image_overlay_plot
; PURPOSE:
;   Make bitmapped overlay for rgb images.
; CALLING SEQUENCE:
;   hogg_image_overlay_plot, xx,yy,naxis1,naxis2,overlay, $
;                            [extras for plot procedure]
; EXAMPLES:
; INPUTS:
;   xx,yy      - points to plot
;   naxis1     - width in pixels to make overlay
;   naxis2     - height in pixels
;   overlay    - overlay to be added to (fine if not set)
;   [extras]   - plotting inputs, just like for "plot" procedure
; OPTIONAL INPUTS:
;   factor     - integer factor to use for antialiasing; default 2;
;                set to 1 for no antialiasing
; OUTPUT:
;   overlay    - overlay with plot material added
; COMMENTS:
;   - Never makes axes!
; BUGS:
;   - Relies on horrifying UNIX bitmapping code.
;   - Makes insecure intermediate PS file.
; DEPENDENCIES:
;   pstopnm etc
; REVISION HISTORY:
;   2004-02-28  written - Hogg
;-
pro hogg_image_overlay_plot, xx,yy,naxis1,naxis2,overlay,factor=factor, $
                             xstyle=xstyle,xrange=xrange, $
                             ystyle=ystyle,yrange=yrange, $
                             _EXTRA=KeywordsForPlot
prefix= 'tmp_hogg_image_overlay_plot'
naxis1= round(naxis1)
naxis2= round(naxis2)
bangp= !P
bangx= !X
bangy= !Y
set_plot, 'PS'
dpi= floor((naxis1/7.5) > (naxis2/10.0))
xsize= double(naxis1)/double(dpi)
ysize= double(naxis2)/double(dpi)
device, filename=prefix+'.ps',xsize=xsize,ysize=ysize,/inches
!P.MULTI= [0,1,1]
!P.POSITION= [0.,0.,1.,1.]
!X.MARGIN= [0,0]
!X.OMARGIN= [0,0]
!Y.MARGIN= !X.MARGIN
!Y.OMARGIN= !Y.OMARGIN
nw_overlay_range, naxis1,naxis2,xrange,yrange
xstyle= 5
ystyle= 5
plot, xx,yy, $
  xstyle=xstyle,xrange=xrange, $
  ystyle=ystyle,yrange=yrange, $
  _EXTRA=KeywordsForPlot
device, /close
!P= bangp
!X= bangx
!Y= bangy
overlay1= 1.-hogg_image_overlay(prefix+'.ps',naxis1,naxis2,factor=factor)
if keyword_set(overlay) then overlay= overlay+overlay1 $
  else overlay= overlay1
return
end
