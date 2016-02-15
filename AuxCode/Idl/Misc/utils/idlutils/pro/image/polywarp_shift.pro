;+
; NAME:
;   polywarp_shift
; PURPOSE:
;   wrapper on polywarp to do a simple shift
; CALLING SEQUENCE
;   newim= polywarp_rotate( im , shift )
; INPUTS:
;   im       - [nx,ny] original image
;   shift    - [2] x,y shift
; OUTPUTS:
;   newim    - [nx,ny] final image
; COMMENTS:
;   polywarp uses a cubic approximation to the sinc function. 
; REVISION HISTORY:
;   2004-04-15  Written - Blanton (NYU)
;-
function polywarp_shift, im, shift

if NOT keyword_set(degree) then degree = 1
if n_elements(interp) EQ 0 then interp = 2
if n_elements(cubic) EQ 0 then cubic = -0.5

; get input positions
tmp= size(im,/dimensions)
ninx= tmp[0]
niny= tmp[1]
inx= double([0,0,ninx-1,ninx-1])
iny= double([0,niny-1,0,niny-1])

; get reference positions
d2r=!DPI/180.
refx=inx+shift[0]
refy=iny+shift[1]

polywarp, inx,iny,refx,refy,degree,Kx,Ky

imout= poly_2d(im,Kx,Ky,interp,ninx,niny,cubic=cubic,missing=0.)

return, imout

end
