;+
; NAME:
;   dmedsmooth
; PURPOSE:
;   median smooth 
; CALLING SEQUENCE:
;   smooth= dmedsmooth(image, invvar, box=)
; INPUTS:
;   image - [nx, ny] input image
;   invvar - [nx, ny] invverse variance (default all 1.)
;   box - box size for smooth
; OUTPUTS:
;   smooth - [nx, ny] smooth image
; COMMENTS:
;   Doesn't weight by inverse variance: just ignores invvar=0 data.
;   If inverse variance isn't supplied, assumes all data good.
; REVISION HISTORY:
;   11-Jan-2006  Written by Blanton, NYU
;-
;------------------------------------------------------------------------------
function dmedsmooth, image, invvar, box=box

nx=(size(image,/dim))[0]
ny=(size(image,/dim))[1]

; Set source object name
soname=filepath('libdimage.'+idlutils_so_ext(), $
                root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')

smooth=fltarr(nx,ny)
if(NOT keyword_set(invvar)) then $
  invvar=fltarr(nx,ny)+1.
retval=call_external(soname, 'idl_dmedsmooth', float(image), $
                     float(invvar), long(nx), long(ny), long(box), $
                     float(smooth))

return, smooth

end
;------------------------------------------------------------------------------
