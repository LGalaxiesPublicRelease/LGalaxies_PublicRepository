;-----------------------------------------------------------------------
;  Convert cartesian coordinates x[0:2] to spherical coordinates 
;  (phi,theta) on unit sphere. The angles are given in the following
;  ranges:
;    0 <= phi < 360
;    0 <= theta <= 180
;  where theta=0 corresponds to the N pole, and theta=180 is the S pole.
;  If you want to convert the results to RA and DEC, use the following
;  relations:
;    ra = phi
;    dec = 90-theta
pro x_to_angles,x,phi,theta
DRADEG = 180.d0/!dpi

phi=DRADEG*atan(x[1,*],x[0,*])
indx=where(phi lt 0.,count)
if(count gt 0) then phi[indx]=phi[indx]+360.
theta=(90.D)-DRADEG*asin(x[2,*]/sqrt(x[0,*]^2+x[1,*]^2+x[2,*]^2))

end
;-----------------------------------------------------------------------
