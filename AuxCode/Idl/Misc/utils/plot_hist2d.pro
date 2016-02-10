pro plot_hist2d,x,y,$
                 xbin=xbin,ybin=ybin,$
                 xmin=xmin,xmax=xmax,$
                 ymin=ymin,ymax=ymax,$
                 log=log,_extra=extra,hist=hist,$
                 xarray=xarray,yarray=yarray,$
                 point=point,psym=psym,pcol=pcol, $
                 cont=cont, frac=frac, weight=weight


; NAME: plot_hist_2d
;
; PURPOSE: plots a 2D histogram of data pairs with grey-scale, contour + points plot
;
;
; CATEGORY: graphics
;
;
; CALLING SEQUENCE:plot_hist_2d,x,y                 
;                 [xbin=xbin,ybin=ybin,$
;                  xmin=xmin,xmax=xmax,$
;                  ymin=ymin,ymax=ymax,$
;                  log=log,_extra=extra,hist=hist,$
;                  xarray=xarray,yarray=yarray,$
;                  point=point,psym=psym,pcol=pcol,cont=cont, frac=frac]
;
;  
;
; INPUTS: x,y : coordinates of points to be binned
;  
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;          xbin: size of bins in x axis (default 1)
;          ybin: size of bins in y axis (default 1)
;           log: if set then the histogram is logged before plotting
;       _ extra: any keywords accepted by contour
;          hist: returns the histogram
;        xarray: returns the x-axis for the histogram
;        yarray: returns the y-axis for the histogram
;         point: if set then points are plotted when histogram is below
;                point threshold
;          psym: symbol to be used for points
;          pcol: colour to be used for points
;          cont: if set then contour plot rather than grey-scale is used
;          frac: if cont and frac are set then contour levels are plotted
;                where density is at frac fractions of the total number of points
;




  xgood=where(finite(x) eq 1)
  ygood=where(finite(y) eq 1)

  if not keyword_set(xbin) then xbin=1
  if not keyword_set(ybin) then ybin=1
  if not keyword_set(xmin) then xmin=min(x[xgood])
  if not keyword_set(xmax) then xmax=max(x[xgood])
  if not keyword_set(ymin) then ymin=min(y[ygood])
  if not keyword_set(ymax) then yax=max(y[ygood])  
  if not keyword_set(psym) then psym=3
  if not keyword_set(pcol) then pcol=0  
  if not keyword_set(weight) then weight=1
 
  hist=hist_2d(x,y,min1=xmin,min2=ymin,max1=xmax,max2=ymax,bin1=xbin,bin2=ybin)*weight


  sz=size(hist)
  xarray=findgen(sz[1])*xbin+xmin
  yarray=findgen(sz[2])*ybin+ymin

  if keyword_set(cont) then begin 

     if keyword_set(log) then begin
        print,'Plotting Log contour'
        contour,alog10(hist>1),xarray+xbin/2.0,yarray+ybin/2.0,_extra=extra,xst=1,yst=1
     endif else BEGIN
        print,'Plotting non-Log contour'
        contour,hist,xarray+xbin/2.0,yarray+ybin/2.0,_extra=extra
     endelse

  endif else begin
     
     ; a lot of this code taken from icplot
     ; Check if device has scalable pixels
     if ( !d.flags and 1 ) then scale = 1 else scale = 0

     ; Set the number of reserved pens
     if not keyword_set( respen ) then respen = 0

     message, 'this bit doesnt work properly - icplot needs fixing'
     if keyword_set(log) then begin
        icplot,alog10(hist>1),xax=xarray,yax=yarray,/pix,/fill,_extra=extra
     endif else begin
        icplot,hist,xax=xarray,yax=yarray,/pix ;,level=10.^(findgen(10)/3.+1),/fill,_extra=extra
     endelse

  endelse




;----------------------------------------------------------------------

  ix=fix((x-xmin)/xbin)
  iy=fix((y-ymin)/ybin)    
  den=hist(ix,iy)

;----------------------------------------------------------------------
; pick contour with levels set at given fractions of the total number
; of points
;----------------------------------------------------------------------
  
  if keyword_set(frac) and keyword_set(cont) then begin

     print, 'hello 2'
     histsm=smooth(float(hist),3)
     densm=histsm(ix,iy)
     order=sort(densm)
     n=round(n_elements(x)*frac)
     den_frac=interpol(densm[order],findgen(n_elements(x)),n)    
 
     contour,histsm,xarray,yarray,/over,levels=den_frac,c_colo=3,xs=1,ys=1
     
  endif

;----------------------------------------------------------------------
; plot up points in sparse areas
;----------------------------------------------------------------------


if keyword_set(point) then begin     
    use=where(den le point,nuse)
    if nuse gt 0 then oplot,x[use]+xbin/2.0*0.0,y[use],psym=psym,col=pcol
endif

 ;contour,alog10(hist>1),xarray+xbin/2.0,yarray+ybin/2.0,_extra=extra,xst=1,yst=1

end
