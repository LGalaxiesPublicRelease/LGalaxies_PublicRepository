;+
; NAME:
;   hogg_histogram
; PURPOSE:
;   Multi-dimensional histogramming function.
; CALLING SEQUENCE:
;   hist = hogg_histogram(data,range,nbin,weight=weight,err=err)
; INPUTS:
;   data     - [d,M] array of d-dimensional data values
;   range    - [2,d] min and max value for every dimension
;   nbin     - [d] array of numbers of bins in each direction
; OPTIONAL INPUTS:
;   weight   - [M] array of weights for the data points
; OUTPUTS:
;   hist     - [d,P] array of numbers of points (or total weight) in each bin
; OPTIONAL OUTPUTS:
;   err      - [d,P] array of Poisson errors
; BUGS:
;   Slow -- very slow.
;   Not memory-efficient.
;   Doesn't deal well with small bin sizes.
;   Doesn't output grid bin positions.
; COMMENTS:
;   Doesn't divide by bin "volumes".
; REVISION HISTORY:
;   2003-07-17  written - Hogg (NYU and NYC Criminal Court)
;-
function hogg_histogram, data,range,nbin,weight=weight,err=err

; check and count inputs
dd= n_elements(nbin)
mm= n_elements(data)/dd
if (n_elements(range) NE (dd*2)) then begin
    splog, 'ERROR: input dimensions do not match'
    return, -1
endif

; reform for safety -- this step may be slow
xx= reform([data],dd,mm)
rr= reform([range],2,dd)
nb= reform([round(nbin)],dd)

; compute binning function -- this step is a memory hog
minarray= transpose(range[0,*])#(dblarr(mm)+1.0)
maxarray= transpose(range[1,*])#(dblarr(mm)+1.0)
nbinarray= nbin#(dblarr(mm)+1L)
bin= djs_floor(nbinarray*(xx-minarray)/(maxarray-minarray))
good= where(total((bin LT 0) OR (bin GE nbinarray),1) EQ 0,ngood) 

; get binning into one-dimensional form -- this step is an IDL-specific hack
onedbin= lonarr(mm)+bin[dd-1,*]
for jj=dd-2L,0L,-1L do onedbin= onedbin*nbin[jj]+bin[jj,*]

; make output array -- check out this IDL hack
hist= reform(dblarr(round(exp(total(alog(nbin),/double)))),nbin)
err= hist

; populate the array -- this step is slow
for kk=0L,ngood-1L do begin
    ii= good[kk]
    if keyword_set(weight) then ww= weight[ii] else ww= 1D0
    hist[onedbin[ii]]= hist[onedbin[ii]]+ww
    err[onedbin[ii]] = err[onedbin[ii]]+ww^2
endfor
err= sqrt(err)

return, hist
end
