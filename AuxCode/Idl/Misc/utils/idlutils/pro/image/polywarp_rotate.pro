;+
; NAME:
;   polywarp_rotate
; PURPOSE:
;   wrapper on polywarp to do a simple rotation
; CALLING SEQUENCE
;   newim= polywarp_rotate( im , theta,  [center=] )
; INPUTS:
;   im       - [nx,ny] original image
;   theta    - angle of rotation (counter-clockwise is position)
; OPTIONAL INPUTS:
;   center    - [2] x/y center to rotate about 
;                (default [(nx-1)*0.5, (ny-1)*0.5])
; OUTPUTS:
;   newim    - [nx,ny] final image
; COMMENTS:
;   polywarp uses a cubic approximation to the sinc function. 
;   this routine rotates the image about its center. 
; REVISION HISTORY:
;   2004-04-15  Written - Blanton (NYU)
;-
function polywarp_rotate, im, theta, center=center

if NOT keyword_set(degree) then degree = 1
if n_elements(interp) EQ 0 then interp = 2
if n_elements(cubic) EQ 0 then cubic = -0.5
if n_elements(shift) NE 2 then shift=[0.,0.]

; get input positions
tmp= size(im,/dimensions)
ninx= tmp[0]
niny= tmp[1]
inx= double([0,0,ninx-1,ninx-1])
iny= double([0,niny-1,0,niny-1])
if(NOT keyword_set(center) ) then $
  center=[float(ninx-1L)*0.5, float(niny-1L)*0.5]

; get reference positions
d2r=!DPI/180.
refx=cos(d2r*theta)*(inx-center[0])-sin(d2r*theta)*(iny-center[1])+center[0]
refy=sin(d2r*theta)*(inx-center[0])+cos(d2r*theta)*(iny-center[1])+center[1]

polywarp, inx,iny,refx,refy,degree,Kx,Ky

imout= poly_2d(im,Kx,Ky,interp,ninx,niny,cubic=cubic,missing=0.)

return, imout

end
