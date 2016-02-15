;+
; NAME:
;   hogg_iterated_ad2xy
; PURPOSE:
;   Same as ad2xy, but assuming that xy2ad is perfect.
; COMMENTS:
;   Starts with a first guess that everything is at (x,y)=(0,0) and
;     then uses Newton's method to iterate to the correct (x,y)
;     values.
; INPUTS:
;   a,d    - vectors of RA,Dec
;   astr   - astrometric structure
; OPTIONAL INPUTS:
;   tol    - maximum allowed angular difference (deg) between input
;            a,d values and the values obtained by running xy2ad on
;            the output; default to 0.01/3600.0 (0.01 arcsec)
; OUTPUTS:
;   x,y    - vectors of x,y
; REVISION HISTORY:
;   2005-09-04  started - Hogg (NYU)
;-
pro hogg_iterated_ad2xy, aa,dd,astr,xx,yy,tol=tol
if (NOT keyword_set(tol)) then tol= 0.01/3600.0
splog, 'making first guess'
xx= dblarr(n_elements(aa))
yy= xx
xy2ad, xx,yy,astr,aaaa,dddd
maxdist= max(djs_diff_angle(aa,dd,aaaa,dddd))
kk= 0
splog, 'iteration',kk, $
  ': get maximum angular error of',maxdist,' (tol',tol,')'
while (maxdist GT tol) do begin
    splog, 'taking derivatives'
    xy2ad, xx+1D0,yy,astr,dadx,dddx
    dadx= temporary(dadx)-aaaa
    dddx= temporary(dddx)-dddd
    xy2ad, xx,yy+1D0,astr,dady,dddy
    dady= temporary(dady)-aaaa
    dddy= temporary(dddy)-dddd
    splog, 'inverting matrices'
    determ= dadx*dddy-dddx*dady
    dxda= dddy/determ
    dxdd= -dady/determ
    dyda= -dddx/determ
    dydd= dadx/determ
    splog, 'updating x,y'
    xx= temporary(xx)+(aa-aaaa)*dxda+(dd-dddd)*dxdd
    yy= temporary(yy)+(aa-aaaa)*dyda+(dd-dddd)*dydd
    xy2ad, xx,yy,astr,aaaa,dddd
    maxdist= max(djs_diff_angle(aa,dd,aaaa,dddd))
    kk= kk+1
    splog, 'iteration',kk, $
      ': get maximum angular error of',maxdist,' (tol',tol,')'
endwhile
return
end
