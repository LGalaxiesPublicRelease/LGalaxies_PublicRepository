;+
; NAME:
;  xyz_to_angles
; PURPOSE:
;  Convert Cartesian coordinates (x,y,z) to spherical coordinates
;  (r,phi,theta).  The returned angles are in the following ranges:
;    0 <= phi < 360
;    0 <= theta <= 180
;  where theta=0 corresponds to the N pole, and theta=180 is the S pole.
;  Note that RA=phi and DEC=90-theta.
; BUGS:
;  - May have divide by zero issues. But why would anyone run this
;    with r=0??
; REVISION HISTORY:
;  2005-09-10  tweaked - Hogg (NYU)
;-
pro xyz_to_angles,x,y,z,r,phi,theta
DRADEG = 180.d0/!dpi
r = sqrt(x*x+y*y+z*z)
theta = DRADEG * acos(z/r)
phi= DRADEG * atan(y,x)
bad= where(phi LT 0.0,nbad)
if (nbad GT 0) then phi[bad] = phi[bad]+3.6D2
return
end
