;+
; NAME:
;   hogg_weighted_mean_surface
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
;   nearest   - use nearest few points based on this number, rather 
;               than a fixed smoothing (use with caution!)
; OUTPUTS:
;   images    - [nx,ny,4] output images of number of contributing
;               points (image 0), total weight used (image 1), total
;               square weight used (image 2), and weighted mean (image 3)
; COMMENTS:
; BUGS:
; REVISION HISTORY:
;   2003-01-11  written - Hogg
;-
function hogg_weighted_mean_surface, x,y,quantity,weight,xbin,ybin,dx,dy, $
                                     boot_seed=boot_seed, nearest=nearest

; set dimensions, defaults, arrays, etc
npt= n_elements(x)
nx= n_elements(xbin)
ny= n_elements(ybin)
image= dblarr(nx,ny,5)
splog, 'making ',nx,' by ',ny,' surface from ',npt,' data points...'

; bootstrap?
if keyword_set(boot_seed) then begin
    boot_seed= long(boot_seed)
    boot_index= floor(randomu(boot_seed,npt)*double(npt))
    tx= x[boot_index]
    ty= y[boot_index]
    tq= quantity[boot_index]
    tw= weight[boot_index]
endif else begin
    tx= x
    ty= y
    tq= quantity
    tw= weight
endelse

; loop over x,y subsamples
if(NOT keyword_set(nearest)) then begin
    for xi= 0L,nx-1 do for yi= 0L,ny-1 do begin
        isub= where(abs(tx-xbin[xi]) LE abs(0.5*dx) AND $
                    abs(ty-ybin[yi]) LE abs(0.5*dy),nsub)
        image[xi,yi,0]= double(nsub)
        if nsub GT 0 then begin
            image[xi,yi,1]= total(tw[isub],/double)
            image[xi,yi,2]= total((tw[isub])^2,/double)
            image[xi,yi,3]= total(tw[isub]*tq[isub],/double)/image[xi,yi,1]
            image[xi,yi,4]=  $
              total(tw[isub]*(tq[isub]-image[xi,yi,3])^2,/double)/ $
              image[xi,yi,1]
        endif
    endfor
endif else begin
    nsub=long(nearest)
    for xi= 0L,nx-1 do for yi= 0L,ny-1 do begin
        dist2=((tx-xbin[xi])^2+(ty-ybin[yi])^2) 
        isort=sort(dist2)
        isub=isort[0:nsub-1L]
        if(sqrt(dist2[isub[0]]) lt dx) then begin
            image[xi,yi,0]= double(nsub)
            image[xi,yi,1]= total(tw[isub],/double)
            image[xi,yi,2]= total((tw[isub])^2,/double)
            image[xi,yi,3]= total(tw[isub]*tq[isub],/double)/image[xi,yi,1]
            image[xi,yi,4]=  $
              total(tw[isub]*(tq[isub]-image[xi,yi,3])^2,/double)/ $
              image[xi,yi,1]
        endif
    endfor
endelse 

; return
splog, '...done'
return, image
end
