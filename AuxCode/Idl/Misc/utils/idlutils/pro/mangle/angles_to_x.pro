;-----------------------------------------------------------------------
;  Convert spherical coordinates (phi,theta) on unit sphere into
;  Cartesion coordinates x[0:2];  The angles must be in the following 
;  ranges:
;    0 <= phi < 360
;    0 <= theta <= 180
;  where theta=0 corresponds to the N pole, and theta=180 is the S pole.
;  If you want to convert from RA and DEC, pass the following
;  arguments (in degrees):  RA, DEC+90
 
function angles_to_x,phi,theta
   DRADEG = 180.d0/!dpi
 
   stheta = sin(theta / DRADEG)
   x=dblarr(3,n_elements(phi))
   x[0,*] = cos(phi / DRADEG) * stheta
   x[1,*] = sin(phi / DRADEG) * stheta
   x[2,*] = cos(theta / DRADEG)
 
   return,x
end
;-----------------------------------------------------------------------
