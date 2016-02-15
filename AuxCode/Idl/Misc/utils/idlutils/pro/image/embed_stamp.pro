;+
; NAME:
;   embed_stamp
; PURPOSE:
;   Add to one image the parts which are overlapped by a second image 
; CALLING SEQUENCE:
;   embed_stamp,full_image,stamp_image,xlo,ylo
; INPUTS:
;   full_image - [full_nx, full_ny] image to add to
;   stamp_image - [stamp_nx, stamp_ny] image to add by
;   xlo, ylo - position in full_image of lower left corner of lower left
;              pixel of stamp
; OUTPUTS:
;   full_image - resulting image
; REVISION HISTORY:
;   2003-01-20  Written - Blanton
;-
pro embed_stamp, full_image, stamp_image, xlo, ylo

if(n_params() ne 4) then begin
    print,'Syntax - embed_stamp, full_image, stamp_image, xlo, ylo'
    return
endif

; sinc-shift stamp_image the fractional part
ixlo=long(xlo)
iylo=long(ylo)
fxlo=xlo-double(ixlo)
fylo=ylo-double(iylo)

; find overlapping region and add
full_nx=(size(full_image,/dimensions))[0]
full_ny=(size(full_image,/dimensions))[1]
shift_stamp_nx=(size(stamp_image,/dimensions))[0]
shift_stamp_ny=(size(stamp_image,/dimensions))[1]
full_xlo=ixlo > 0
full_xhi=(ixlo+shift_stamp_nx) < full_nx
full_ylo=iylo > 0
full_yhi=(iylo+shift_stamp_ny) < full_ny
shift_stamp_xlo=full_xlo-ixlo
shift_stamp_xhi=full_xhi-ixlo
shift_stamp_ylo=full_ylo-iylo
shift_stamp_yhi=full_yhi-iylo
if(shift_stamp_xlo ge 0 and shift_stamp_ylo ge 0 and $
   shift_stamp_xhi gt shift_stamp_xlo and $
   shift_stamp_yhi gt shift_stamp_ylo and $
   full_xhi gt full_xlo and full_yhi gt full_ylo) then begin
    tmp_image=stamp_image[shift_stamp_xlo:shift_stamp_xhi-1L, $
                          shift_stamp_ylo:shift_stamp_yhi-1L]
    tmp_image=sshift2d(tmp_image,[fxlo,fylo])
    full_image[full_xlo:full_xhi-1L,full_ylo:full_yhi-1L]= $
      full_image[full_xlo:full_xhi-1L,full_ylo:full_yhi-1L]+tmp_image
endif

return
end
