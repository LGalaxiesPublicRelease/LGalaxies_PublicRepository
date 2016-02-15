;+
; NAME:
;   dsmooth
; PURPOSE:
;   smooth with a simple gaussian
; CALLING SEQUENCE:
;   smooth= dsmooth(image, sigma)
; INPUTS:
;   image - [nx, ny] input image
;   sigma - gaussian sigma
; OUTPUTS:
;   smooth - [nx, ny] smooth image
; REVISION HISTORY:
;   11-Jan-2006  Written by Blanton, NYU
;-
;------------------------------------------------------------------------------
function dsmooth, image, sigma

if((size(image))[0] eq 1) then begin
    ny=1
    nx=(size(image,/dim))[0]
endif else begin
    nx=(size(image,/dim))[0]
    ny=(size(image,/dim))[1]
endelse

; Set source object name
soname=filepath('libdimage.'+idlutils_so_ext(), $
                root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')

smooth=fltarr(nx,ny)
retval=call_external(soname, 'idl_dsmooth', float(image), $
                     long(nx), long(ny), float(sigma), float(smooth))

return, smooth

end
;------------------------------------------------------------------------------
