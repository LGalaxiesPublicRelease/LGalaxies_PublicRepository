;+
; NAME:
;  hogg_open_overlay
; PURPOSE:
;  open Z buffer for writing an image overlay plot
; INPUTS:
;  naxis        - [naxis1,naxis2] size of image
; KEYWORDS:
;  noantialias  - by default this uses very slow and memory-hogging
;                 resampling for anti-aliasing
; BUGS:
;  - Header not complete.
;  - Line thicknesses and other crap hard-coded.
;  - !P,!X,!Y variable over-written and not saved; keep in common block.
;-
pro hogg_open_overlay, naxis,noantialias=noantialias
set_plot,'Z'
if keyword_set(noantialias) then factor=1 else factor=2
resolution= naxis*factor
device, set_resolution=resolution,z_buffer=0
hogg_plot_defaults
!P.THICK= 12.0*naxis[0]/2048.0
!P.CHARSIZE= 2.0*!P.THICK
!P.CHARTHICK= 2.0*!P.THICK
!X.OMARGIN=[0,0]
!Y.OMARGIN=[0,0]
loadct, 0
!P.BACKGROUND= 255
!P.COLOR= 0
plot, [0,1],/nodata, $
  xrange=[0,naxis[0]]-0.5,xstyle=5, $
  yrange=[0,naxis[1]]-0.5,ystyle=5
return
end
