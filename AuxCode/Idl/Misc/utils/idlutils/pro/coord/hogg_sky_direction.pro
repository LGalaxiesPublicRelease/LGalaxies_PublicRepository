;+
; NAME:
;   hogg_sky_direction
; PURPOSE:
;   Return the direction (Delta-alpha,Delta-delta) of one RA, Dec
;   relative to a reference RA, Dec.
; CALLING SEQUENCE:
;   hogg_sky_direction, aref,dref,ap,dp,Da,Dd
; INPUT:
;   aref,dref  - reference RA,Dec
;   ap,dp      - RA,Dec of the point to which the arrow points
; OUTPUT:
;   Da         - Delta-alpha (small change in RA*cos(Dec)) for arrow
;   Dd         - Delta-delta (small change in Dec) for arrow
; COMMENTS:
;   - Da,Dd comprise a vector of arbitrary length; it is the user's
;     responsibility to re-scale it.
;   - Da is not a change in RA, but a change in RA*cos(Dec).
; BUGS:
;   - Inputs must all be scalars or else hell breaks loose.
; REVISION HISTORY:
;   2005-09-10  started - Hogg (NYU)
;-
pro hogg_sky_direction, aref,dref,ap,dp,Da,Dd

; make unit vector to reference point
angles_to_xyz, 1D0,aref,9D1-dref,xref,yref,zref
ref= [xref,yref,zref]

; make unit vector to object point and then perpendicular-project it!
angles_to_xyz, 1D0,ap,9D1-dp,xp,yp,zp
pp= [xp,yp,zp]
dotprod= (transpose(ref)#pp)[0]
pp= pp-dotprod*ref ; pp is now perpendicular to ref
pp= pp/sqrt(total(pp^2))

; perpendicular-project zhat to get the increasing-delta direction
zhat= [0D0,0D0,1D0]
dotprod= (transpose(ref)#zhat)[0]
deltahat= zhat-dotprod*ref
deltahat= deltahat/sqrt(total(deltahat^2)) ; delta direction

; make the increasing-alpha direction with a cross product
alphahat= [deltahat[1]*ref[2]-deltahat[2]*ref[1], $
           deltahat[2]*ref[0]-deltahat[0]*ref[2], $
           deltahat[0]*ref[1]-deltahat[1]*ref[0]] ; alpha direction

; project and return
Da= (transpose(alphahat)#pp)[0]
Dd= (transpose(deltahat)#pp)[0]
return
end
