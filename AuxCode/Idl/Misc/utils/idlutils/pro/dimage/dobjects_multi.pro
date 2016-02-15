;+
; NAME:
;   dobjects_multi
; PURPOSE:
;   detect objects in a set of images of different bands
; CALLING SEQUENCE:
;   dobjects_multi, images, objects= [, dpsf=, plim=]
; INPUTS:
;   images - [nx, ny, nim] original image
; OPTIONAL INPUTS:
;   dpsf - smoothing of PSF for detection (defaults to sigma=1 pixel)
;   plim - limiting significance in sky sigma (defaults to 10 sig )
; OUTPUTS:
;   objects - [nx, ny] which object each pixel belongs to (-1 if none)
; COMMENTS:
;   Calculates noise in image from image itself.
;   Any detected pixel in any band counts as a detection.
; REVISION HISTORY:
;   11-Jan-2006  Written by Blanton, NYU
;-
;------------------------------------------------------------------------------
pro dobjects_multi, images, objects=objects, dpsf=dpsf, plim=plim

nx=(size(images,/dim))[0]
ny=(size(images,/dim))[1]
if((size(images))[0] eq 2) then $
  nim=1 $
else $
  nim=(size(images,/dim))[2]

if(NOT keyword_set(dpsf)) then dpsf=1.
if(NOT keyword_set(plim)) then plim=10.

; Set source object name
soname=filepath('libdimage.'+idlutils_so_ext(), $
                root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')

nchild=0L
objects=lonarr(nx,ny)
retval=call_external(soname, 'idl_dobjects_multi', $
                     float(images), $
                     long(nx), long(ny), long(nim), $
                     float(dpsf), $
                     float(plim), $
                     long(objects))

end
