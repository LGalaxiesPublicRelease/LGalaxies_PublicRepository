;+
; NAME:
;  angles_to_xyz
; PURPOSE:
;  Convert spherical coordinates (r,phi,theta) into Cartesion coordinates
;  (x,y,z).  The angles must be in the following ranges:
;    0 <= phi < 360
;    0 <= theta <= 180
;  where theta=0 corresponds to the N pole, and theta=180 is the S pole.
;  If you want to convert from RA and DEC, pass the following
;  arguments (in degrees):  RA, 90-DEC
; REVISION HISTORY:
;  2005-09-10  tweaked - Hogg (NYU)
;-
pro angles_to_xyz,r,phi,theta,x,y,z
DRADEG = 180.d0/!dpi
x = r * cos(phi / DRADEG) * sin(theta / DRADEG)
y = r * sin(phi / DRADEG) * sin(theta / DRADEG)
z = r * cos(theta / DRADEG)
return
end
