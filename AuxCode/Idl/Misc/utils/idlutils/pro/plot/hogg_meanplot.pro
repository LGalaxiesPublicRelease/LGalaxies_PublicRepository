;+
; NAME:
;   hogg_meanplot
; PURPOSE:
;   plot sliding mean of one quantity vs two others
; COMMENTS:
;   Doesn't overplot -- only plots.  This is because IDL blots
;     out any other information on the plot anyway.
; INPUTS:
;   x,y         - data values
;   z           - quantity to average
; OPTIONAL INPUTS:
;   weight      - weighting for data points; default unity
;   xrange      - x range; default to the number > minnum area
;   yrange      - y range; same default
;   dxbin       - size of boxes in x-dir; default to a function of
;                 first and second moments 
;   dybin       - size of boxes in y-dir; same default
;   levels      - contour levels; default to a function of image range
;   minnum      - minimum number of points in a sliding box to plot;
;                 default 100
;   c_colors    - fill colors
;   c_thick     - contour line thicknesses
;   axis_char_scale - scale to apply to axis labels
;   nearest     - take mean over nearest few points based on this number
;                 (instead of a fixed smoothing scale)
; KEYWORDS:
;   nofill      - don't fill the contours with color
;   noperimeter - don't plot contour at minnum
;   nobox       - don't plot box
;   nolines     - don't plot lines between contours
;   nodata      - don't plot anything other than axes
;   maskonly    - plot only the minnum mask
;   nocontourlabels  - don't label the contours
;   input_mean  - don't recalculate means, use input in bin_mean
;   overplot    - don't remake axes, just plot over
; OPTIONAL INPUT/OUTPUTS:
;   bin_mean    - calculated mean value (or input mean value if 
;                 /input_mean set) in each bin
;   bin_scatter - calculated scatter around the mean in each bin
;   bin_weight  - total weight in each bin
;   bin_number  - total number of data points in each bin
; COMMENTS:
;   "nearest" input should be used with caution --- this is generally
;     a BAD way to smooth!
; REVISION HISTORY:
;   2003-01-08  written - Hogg
;-
pro hogg_meanplot, x,y,z,weight=weight, $
                   xrange=xrange,yrange=yrange,dxbin=dxbin,dybin=dybin, $
                   levels=levels,c_colors=c_colors,c_thick=c_thick, $
                   minnum=minnum, nodata=nodata, $
                   noperimeter=noperimeter,nobox=nobox,nolines=nolines, $
                   maskonly=maskonly, xbin=xbin, ybin=ybin, bin_mean=bin_mean, $
                   bin_scatter=bin_scatter, bin_weight=bin_weight, $
                   bin_number=bin_number, input_mean=input_mean, $
                   axis_char_scale=axis_char_scale, $
                   nocontourlabels=nocontourlabels, overplot=overplot, $
                   nofill=nofill, nearest=nearest

if(NOT keyword_set(nofill)) then cell_fill=1L
if(NOT keyword_set(minnum)) then minnum=1L
if(NOT keyword_set(axis_char_scale)) then axis_char_scale=1.75

; take moments
ndata= n_elements(x)
if not keyword_set(weight) then weight= fltarr(ndata)+1.0
weightsum= total(weight,/double)
qmean= total(weight*z,/double)/weightsum
qvar= total(weight*(z-qmean)^2,/double)/weightsum
xmean= total(weight*x,/double)/weightsum
xvar= total(weight*(x-xmean)^2,/double)/weightsum
ymean= total(weight*y,/double)/weightsum
yvar= total(weight*(y-ymean)^2,/double)/weightsum

; set defaults
if not keyword_set(dxbin) then dxbin= sqrt(xvar)/3.0
if not keyword_set(dybin) then dybin= sqrt(yvar)/3.0
if not keyword_set(minnum) then minnum= 100

; deal with negatives harshly
dxbin= abs(dxbin)
dybin= abs(dybin)

; make (overlapping) bins
factor= 3.0                     ; number of bin centers per bin width
factor2= 8.0                    ; number of sigmas to cover
if not keyword_set(xrange) then begin ; base the bins on moments if no xrange
    nxbin= ceil(factor*2.0*factor2*sqrt(xvar)/dxbin)+1
    xbin= xmean+(dxbin/factor)*(dindgen(nxbin)-0.5*double(nxbin-1))
endif else begin                ; base the bins on xrange if possible
    nxbin= ceil(factor*abs((xrange[1]-xrange[0])/dxbin))+1
    sign= (xrange[1]-xrange[0])/abs(xrange[1]-xrange[0])
    xbin= xrange[0]+sign*(dxbin/factor)*dindgen(nxbin)
endelse
if not keyword_set(yrange) then begin
    nybin= ceil(factor*2.0*factor2*sqrt(yvar)/dybin)+1
    ybin= ymean+(dybin/factor)*(dindgen(nybin)-0.5*double(nybin-1))
endif else begin
    nybin= ceil(factor*abs((yrange[1]-yrange[0])/dybin))+1
    sign= (yrange[1]-yrange[0])/abs(yrange[1]-yrange[0])
    ybin= yrange[0]+sign*(dybin/factor)*dindgen(nybin)
endelse

; make mean image
if (NOT keyword_set(input_mean)) then begin
    image= hogg_weighted_mean_surface(x,y,z,weight,xbin,ybin,dxbin,dybin, $
                                      nearest=nearest)
    bin_number= image[*,*,0]
    bin_weight= image[*,*,1]
    bin_mean= image[*,*,3]
    bin_scatter= image[*,*,4]
endif

; check values and set contour levels
factor= 10.0
good= where(bin_number GT minnum,count_good)
if(count_good gt 0) then begin
    limits_indx=good
endif else begin
    splog,'WARNING: no bins above minnum'
    limits_indx=lindgen(n_elements(bin_number))
endelse
temp= minmax(bin_mean[limits_indx])
if not keyword_set(levels) then $
  levels= temp[0]+(temp[1]-temp[0])/factor*dindgen(ceil(factor+1.0))
nlevels= n_elements(levels)
if (min(levels) EQ max(levels)) then begin
    levels= levels[0]
    nlevels= 1
endif

; set x and y ranges
if not keyword_set(xrange) then begin
    xbin_image= xbin#(dblarr(n_elements(ybin))+1)
    xrange= minmax(xbin_image[limits_indx])
endif
if not keyword_set(yrange) then begin
    ybin_image= (dblarr(n_elements(xbin))+1)#ybin
    yrange= minmax(ybin_image[limits_indx])
endif

; set plot range
nticks=6/axis_char_scale
!X.TICKINTERVAL= hogg_interval(xrange,nticks=nticks)
!Y.TICKINTERVAL= hogg_interval(yrange,nticks=nticks)

; make contour plot
loadct,0,/silent
if NOT keyword_set(c_colors) then begin
    c_colors= 255.0-128.0*(dindgen(nlevels))/double(nlevels)
endif
contour, bin_mean,xbin,ybin,levels=levels,cell_fill=cell_fill, $
  c_colors=c_colors, overplot=overplot, $
  xstyle=1,xrange=xrange,ystyle=1,yrange=yrange, nodata=nodata
if NOT keyword_set(nolines) then begin
    if NOT keyword_set(nocontourlabels) then begin
        c_labels= lonarr(n_elements(levels))+1L
    endif else begin
        c_labels= lonarr(n_elements(levels))
    endelse
    contour, bin_mean,xbin,ybin,levels=levels,/overplot, $
      c_labels=c_labels,c_charthick=!P.CHARTHICK, $
      nodata=nodata,c_thick=c_thick
endif

; show sliding box
if not keyword_set(nobox) then begin
    box_x= 0.5*[-1,-1,1,1,-1]
    box_y= 0.5*[-1,1,1,-1,-1]
    box_xc= mean(!X.CRANGE)
    box_yc= mean(!Y.CRANGE)
    oplot, dxbin*box_x+box_xc,dybin*box_y+box_yc,psym=0,thick=0.5*!P.THICK
endif

; plot number perimeter
numlevels=[0,minnum]
if not keyword_set(noperimeter) and not keyword_set(nodata) then $
  contour, bin_number,xbin,ybin,levels=numlevels, $
  /overplot,cell_fill=maskonly, c_thick=3

end
