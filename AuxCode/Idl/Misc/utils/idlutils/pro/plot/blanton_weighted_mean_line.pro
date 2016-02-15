;+
; NAME:
;   blanton_weighted_mean_surface
; PURPOSE:
;   make an image of the weighted mean for a 2-d set of measurements
; CALLING SEQUENCE:
;   images = weighted_mean_surface(x,y,weight,quantity,xbin,ybin,dx,dy)
; INPUTS:
;   x,y       - [N] values of measurements
;   quantity  - [N] measurements themselves
;   weight    - [N] weights for measurements
;   xbin,ybin - [nx],[ny] vectors of coordinates of image pixel centers
;   dx,dy     - size of sliding box in which means are taken
; OPTIONAL INPUTS:
;   boot_seed - if set, use as seed for a bootstrap resampling trial
; OUTPUTS:
;   images    - [nx,ny,4] output images of number of contributing
;               points (image 0), total weight used (image 1), total
;               square weight used (image 2), and weighted mean (image 3)
; COMMENTS:
; BUGS:
; REVISION HISTORY:
;   2003-01-11  written - Hogg
;-
function blanton_weighted_mean_line, x,quantity,weight,xbin,dx, $
                                     boot_seed=boot_seed

; set dimensions, defaults, arrays, etc
npt= n_elements(x)
nx= n_elements(xbin)
image= dblarr(nx,5)
splog, 'making ',nx,' surface from ',npt,' data points...'

; bootstrap?
if keyword_set(boot_seed) then begin
    boot_seed= long(boot_seed)
    boot_index= floor(randomu(boot_seed,npt)*double(npt))
    tx= x[boot_index]
    tq= quantity[boot_index]
    tw= weight[boot_index]
endif else begin
    tx= x
    tq= quantity
    tw= weight
endelse

; loop over x,y subsamples
for xi= 0L,nx-1 do begin
    isub= where(abs(tx-xbin[xi]) LE abs(0.5*dx),nsub)
    image[xi,0]= double(nsub)
    if nsub GT 0 then begin
        image[xi,1]= total(tw[isub],/double)
        image[xi,2]= total((tw[isub])^2,/double)
        image[xi,3]= total(tw[isub]*tq[isub],/double)/image[xi,1]
        image[xi,4]= total(tw[isub]*(tq[isub]-image[xi,3])^2,/double)/ $
          image[xi,1]
    endif
endfor

; return
splog, '...done'
return, image
end
