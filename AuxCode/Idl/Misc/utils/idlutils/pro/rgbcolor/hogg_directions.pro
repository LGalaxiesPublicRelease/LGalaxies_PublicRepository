;+
; NAME:
;   hogg_directions
; PURPOSE:
;   Overplot directions to other sky locations on a picture of the
;   sky.
; INPUTS:
;   aref,dref  - reference RA,Dec
;   aa,dd      - RA,Dec lists for things to point to
;   name       - list of names of those things
;   hdr        - FITS header with relevant astrometry structure
; OPTIONAL INPUTS:
;   length     - length of arrows in degrees; default 1/60
; OUTPUTS:
;   [overlay on currently open plot]
; BUGS:
;   - Not tested.
;   - Doesn't work with GSSS headers.
;   - Hacks from nw_ad_grid duplicated here, stupidly.
; REVISION HISTORY:
;   2005-09-11  started - Hogg (NYU)
;-
pro hogg_directions, aref,dref,aa,dd,name,hdr, $
                     length=length,linethick=linethick,charsize=charsize
RADEG= 1.8D2/!DPI
narrow= n_elements(aa)
nx= sxpar(hdr,'NAXIS1')
if (NOT keyword_set(length)) then length= 1D0/6D1
IF (NOT keyword_set(linethick)) THEN linethick= 4
lthick= 2000.*linethick/NX ; empirical HACK
IF (NOT keyword_set(charsize)) THEN charsize= lthick/2.0
charthick= charsize*3.0
extast, hdr,astr
ad2xy, aref,dref,astr,x0,y0
for ii=0L,narrow-1L do begin
    hogg_sky_direction, aref,dref,aa[ii],dd[ii],Deltaa,Deltad
    ad2xy, aref+length*Deltaa/cos(dref/RADEG),dref+length*Deltad,astr,x1,y1
    if (x1 GT x0) then align=0 else align=1
    hogg_arrow, x0,y0,x1,y1,/data,head_angle=10.0, $
      thick=lthick
    xyouts, x1,y1,name[ii],align=align, $
      charsize=charsize,charthick=charthick
endfor
return
end
